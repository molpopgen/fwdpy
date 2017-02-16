/* Write output from fwdpy::selected_mut_tracker
 * to sqlite data bases.
 *
 * This files represents first attempts at learning the
 * sqlite3 C API in C++.
 *
 * Goal:
 * We want to write either to different databases/thread
 * or to 1/db from all threads.  The latter requires
 * synchronization across threads, which is done via std::mutex
 * and std::lock_guard.
 *
 * Solution:
 * We have a base class that is capable of writing to different
 * db/thread.  A derived class redefines the relevant member functions
 * and uses the locks in order to write to 1 db from many threads.
 *
 * Note:
 * The class design itself is likely sub-optimal in a few respects,
 * but it works.
 *
 * Lessons learned:
 * 1. There are a lot of places where error-checking of sqlite3
 * return values occurs.  Given that we want to throw an exception
 * when return values indicate a problem, it is better to encapsulate
 * the operations into a class rather than rely on a large set of functions.
 *
 * 2. This file is part of fwdpy, which is a shared object in a Python module.
 * The C++ standard does not cover the behavior of dynamic libraries, meaning
 * that it
 * is conceivable that a global mutex variable in this file could affect ALL
 * instances
 * of fwdpy being run on the same machine.  Thus, we pass in shared_ptr<mutex>
 * to all threads
 * from the calling environment (which in this case is Python via Cython).
 */
#include <future>
#include <memory>
#include <mutex>
#include <functional>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cstdio>
#include <sqlite3.h>
#include "sampler_selected_mut_tracker.hpp"

using namespace std;

namespace
{
    int
    apply_sql_pragma(sqlite3 *db, char *error_message)
    {
        // from http://blog.quibb.org/2010/08/fast-bulk-inserts-into-sqlite/
        int rc = sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL,
                              &error_message);

        if (rc != SQLITE_OK)
            {
                return rc;
            }
        rc = sqlite3_exec(db, "PRAGMA count_changes=OFF", NULL, NULL,
                          &error_message);
        if (rc != SQLITE_OK)
            {
                return rc;
            }
        rc = sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL,
                          &error_message);
        if (rc != SQLITE_OK)
            {
                return rc;
            }
        rc = sqlite3_exec(db, "PRAGMA cache_size(10000)", NULL, NULL,
                          &error_message);
        if (rc != SQLITE_OK)
            {
                return rc;
            }
        rc = sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL,
                          &error_message);
        return rc;
    }

    int
    execute_sql_statement(sqlite3 *db, const string &sql, char *error_message)
    {
        return sqlite3_exec(db, sql.c_str(), NULL, NULL, &error_message);
    }

    static int
    sql_callback(void *data, int argc, char **argv, char **azColName)
    {
        auto dp = reinterpret_cast<sqlite3_stmt *>(data);
        int idx = 1;
        sqlite3_bind_int(dp, idx, atoi(argv[idx - 1])); // generation
        idx++;
        sqlite3_bind_int(dp, idx, atoi(argv[idx - 1]));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx - 1], NULL));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx - 1], NULL));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx - 1], NULL)); // freq
        int rc = sqlite3_step(dp);
        if (rc != SQLITE_DONE)
            {
                throw runtime_error("sql_callback failed: "
                                    + to_string(__LINE__));
            }
        sqlite3_reset(dp);
        return 0;
    };

    static int
    sql_callback_onedb(void *data, int argc, char **argv, char **azColName)
    {
        auto dp = reinterpret_cast<sqlite3_stmt *>(data);
        int idx = 1;
        sqlite3_bind_int(dp, idx, atoi(argv[idx - 1])); // rep
        idx++;
        sqlite3_bind_int(dp, idx, atoi(argv[idx - 1])); // generation
        idx++;
        sqlite3_bind_int(dp, idx, atoi(argv[idx - 1]));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx - 1], NULL));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx - 1], NULL));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx - 1], NULL)); // freq
        int rc = sqlite3_step(dp);
        if (rc != SQLITE_DONE)
            {
                throw runtime_error("sql_callback_onedb failed: "
                                    + to_string(__LINE__));
            }
        sqlite3_reset(dp);
        return 0;
    };

    struct trajSQL
    {
        sqlite3 *db, *memdb;
        sqlite3_stmt *stmt, *memdb_stmt;
        const fwdpy::trajFilter *tf;
        char *error_message;
        unsigned threshold;

        using callback_fxn_t = int (*)(void *, int, char **, char **);
        explicit trajSQL(sqlite3 *db_, const fwdpy::trajFilter *tf_,
                         const string &dbname_, const unsigned threshold_,
                         const bool append_);
        virtual ~trajSQL() noexcept;
        void handle_return_values(int rc, int line_num);
        virtual void create_table(sqlite3 *db_); // needs overloading -- done
        virtual void create_index(sqlite3 *db_); // needs overloading -- done
        virtual void prepare_statements();       // needs overloading -- Done
        virtual void copy_from_mem_db(callback_fxn_t c);
        virtual unsigned
        apply_prepared_statement(sqlite3 *db_, sqlite3_stmt *stmt_,
                                 const unsigned origin, const double pos,
                                 const double esize,
                                 const vector<pair<unsigned, double>> &freqs);
        virtual void
        call_operator_details(const fwdpy::selected_mut_tracker::final_t
                                  &data); // needs overloading -- Done
        void operator()(const fwdpy::selected_mut_tracker::final_t &data);
    };

    trajSQL::trajSQL(sqlite3 *db_, const fwdpy::trajFilter *tf_,
                     const string &dbname_, const unsigned threshold_,
                     const bool append_)
        : db(nullptr), memdb(nullptr), stmt(nullptr), memdb_stmt(nullptr),
          tf(tf_), error_message(nullptr), threshold(threshold_)
    {
        if (!append_)
            {
                remove(dbname_.c_str());
            }
        if (db_ != nullptr)
            {
                db = db_;
            }
        else
            {
                int rc = sqlite3_open(dbname_.c_str(), &db);
                handle_return_values(rc, __LINE__);
                create_table(db);
                create_index(db);
                rc = apply_sql_pragma(db, error_message);
                handle_return_values(rc, __LINE__);
            }
    }

    trajSQL::~trajSQL() noexcept
    {
        if (db != nullptr)
            {
                sqlite3_close(db);
            }
        if (memdb != nullptr)
            {
                sqlite3_close(memdb);
            }
        if (stmt != nullptr)
            {
                sqlite3_finalize(stmt);
            }
        if (memdb_stmt != nullptr)
            {
                sqlite3_finalize(memdb_stmt);
            }
        if (error_message != nullptr)
            {
                sqlite3_free(error_message);
            }
    }

    void
    trajSQL::handle_return_values(int rc, int line_num)
    {
        if (rc != SQLITE_OK)
            {
				if(error_message==nullptr)
				{
					throw runtime_error("error encounted without message: " + to_string(rc) + " " + to_string(line_num));
				}
                string message(error_message);
                message += " ";
                message += to_string(line_num);
                throw runtime_error(message);
            }
    }

    void
    trajSQL::create_table(sqlite3 *db_)
    {
        string sql("CREATE TABLE IF NOT EXISTS freqs(");
        sql += "generation int NOT NULL, origin int NOT NULL, pos real NOT ";
        sql += "NULL, esize real NOT NULL, freq real NOT NULL);";
        int rc = execute_sql_statement(db_, sql, this->error_message);
        handle_return_values(rc, __LINE__);
    }

    void
    trajSQL::create_index(sqlite3 *db_)
    {
        string sql = "create index if not exists gen on freqs (generation);";
        int rc = execute_sql_statement(db_, sql, this->error_message);
        handle_return_values(rc, __LINE__);
    }

    void
    trajSQL::prepare_statements()
    {
        string prepped_statement
            = "insert into freqs values (?1,?2,?3,?4,?5);";
        sqlite3_prepare_v2(db, prepped_statement.c_str(),
                           prepped_statement.size(), &stmt, NULL);
        if (stmt == NULL)
            {
                throw runtime_error(
                    "could not prepare statement for on-disk database.");
            }
        if (memdb != nullptr)
            {
                sqlite3_prepare_v2(memdb, prepped_statement.c_str(),
                                   prepped_statement.size(), &memdb_stmt,
                                   NULL);
                if (memdb_stmt == NULL)
                    {
                        throw runtime_error("could not prepare statement for "
                                            "in-memory database.");
                    }
            }
    }

    unsigned
    trajSQL::apply_prepared_statement(
        sqlite3 *db_, sqlite3_stmt *stmt_, const unsigned origin,
        const double pos, const double esize,
        const vector<pair<unsigned, double>> &freqs)
    {
        if (stmt_ == NULL)
            {
                throw std::runtime_error("NULL statement encountered");
            }
        for (auto &&fi : freqs)
            {
                int idx = 1;
                sqlite3_bind_int(stmt_, idx++,
                                 static_cast<int>(fi.first)); // generation
                sqlite3_bind_int(stmt_, idx++, static_cast<int>(origin));
                sqlite3_bind_double(stmt_, idx++, pos);
                sqlite3_bind_double(stmt_, idx++, esize);
                sqlite3_bind_double(stmt_, idx++, fi.second); // freq
                int rc = sqlite3_step(stmt_);
                if (rc != SQLITE_DONE)
                    {
                        ostringstream message;
                        message << "Prepared statement could not be executed: "
                                << this_thread::get_id() << ' ' << __LINE__
                                << ' ' << idx << ' ' << rc << '\n';
                        throw std::runtime_error(message.str().c_str());
                    }
                sqlite3_reset(stmt_);
            }
        return static_cast<unsigned>(freqs.size());
    }

    void
    trajSQL::copy_from_mem_db(callback_fxn_t c)
    {
        string sql = "select * from freqs";
        int rc = sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL,
                              &error_message);
        handle_return_values(rc, __LINE__);
        rc = sqlite3_exec(memdb, sql.c_str(), c, (void *)stmt, &error_message);
        handle_return_values(rc, __LINE__);
        rc = sqlite3_exec(db, "COMMIT TRANSACTION;", NULL, NULL,
                          &error_message);
        handle_return_values(rc, __LINE__);
        sql = "delete from freqs;";
        rc = sqlite3_exec(memdb, sql.c_str(), NULL, NULL, &error_message);
        handle_return_values(rc, __LINE__);
    }
    void
    trajSQL::call_operator_details(
        const fwdpy::selected_mut_tracker::final_t &data)
    {
        // unsigned nrecords_passed = 0;

        sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &error_message);
        for (auto &&outer : data)
            {
                if (this->tf->apply_origin_filter(outer.first))
                    {
                        for (auto &&inner : outer.second)
                            {
                                if (this->tf->apply_pos_esize_filter(
                                        inner.first)
                                    && (this->tf->apply_freq_filter(
                                           inner.second)))
                                    {
                                        apply_prepared_statement(
                                            db, stmt, outer.first,
                                            inner.first.first,
                                            inner.first.second, inner.second);
                                    }
                            }
                    }
            }
        sqlite3_exec(db, "COMMIT TRANSACTION;", NULL, NULL, &error_message);
    }

    void
    trajSQL::operator()(const fwdpy::selected_mut_tracker::final_t &data)
    {
        call_operator_details(data);
    }

    class trajSQLonedb : public trajSQL
    {
      private:
        shared_ptr<mutex> dblock;
        sqlite3 *init_db_mutex(const string &dbname, const bool append,
                               shared_ptr<mutex> dblock_copy);

      public:
        const unsigned label;
        using base = trajSQL;
        explicit trajSQLonedb(const shared_ptr<mutex> &dblock_,
                              const fwdpy::trajFilter *tf,
                              const string &dbname_, const unsigned threshold_,
                              const unsigned label_, const bool append_);
        void prepare_statements() final;
        unsigned apply_prepared_statement(
            sqlite3 *db_, sqlite3_stmt *stmt_, const unsigned origin,
            const double pos, const double esize,
            const vector<pair<unsigned, double>> &freqs) final;
        void call_operator_details(
            const fwdpy::selected_mut_tracker::final_t &data) final;
        void create_table(sqlite3 *db) final;
        void create_index(sqlite3 *db) final;
    };

    trajSQLonedb::trajSQLonedb(const shared_ptr<mutex> &dblock_,
                               const fwdpy::trajFilter *tf,
                               const string &dbname_,
                               const unsigned threshold_,
                               const unsigned label_, const bool append_)
        : base(init_db_mutex(dbname_, append_, dblock_), tf, dbname_,
               threshold_, true),
          dblock(dblock_), label(label_)
    {
        lock_guard<mutex> lock(*(this->dblock));
        create_table(db);
        create_index(db);
        int rc = apply_sql_pragma(db, error_message);
        handle_return_values(rc, __LINE__);
        rc = sqlite3_open(":memory:", &memdb);
        handle_return_values(rc, __LINE__);
        create_table(memdb);
        rc = apply_sql_pragma(memdb, error_message);
        handle_return_values(rc, __LINE__);
    }

    sqlite3 *
    trajSQLonedb::init_db_mutex(const string &dbname, const bool append,
                                shared_ptr<mutex> dblock_copy)
    {
        if (!append)
            {
                remove(dbname.c_str());
            }
        sqlite3 *rv;
        lock_guard<mutex> lock(*dblock_copy);
        int rc = sqlite3_open(dbname.c_str(), &rv);
        if (rc != SQLITE_OK)
            {
                throw std::runtime_error("Could not open " + dbname);
            }
        return rv;
    }

    void
    trajSQLonedb::create_table(sqlite3 *db_)
    {
        string sql("CREATE TABLE IF NOT EXISTS freqs(");
        sql += "rep int not null, generation int NOT NULL, origin int NOT "
               "NULL, pos real NOT ";
        sql += "NULL, esize real NOT NULL, freq real NOT NULL);";
        int rc = execute_sql_statement(db_, sql, this->error_message);
        handle_return_values(rc, __LINE__);
    }

    void
    trajSQLonedb::create_index(sqlite3 *db_)
    {
        string sql
            = "create index if not exists rep_gen on freqs (rep,generation);";
        int rc = execute_sql_statement(db_, sql, this->error_message);
        handle_return_values(rc, __LINE__);
    }

    void
    trajSQLonedb::prepare_statements()
    {
        string prepped_statement
            = "insert into freqs values (?1,?2,?3,?4,?5,?6);";
        sqlite3_prepare_v2(db, prepped_statement.c_str(),
                           prepped_statement.size(), &stmt, NULL);
        if (stmt == NULL)
            {
                throw runtime_error(
                    "could not prepare statement for on-disk database: "
                    + to_string(__LINE__));
            }
        sqlite3_prepare_v2(memdb, prepped_statement.c_str(),
                           prepped_statement.size(), &memdb_stmt, NULL);
        if (memdb_stmt == NULL)
            {
                throw runtime_error(
                    "could not prepare statement for in-memory database: "
                    + to_string(__LINE__));
            }
    }
    unsigned
    trajSQLonedb::apply_prepared_statement(
        sqlite3 *db_, sqlite3_stmt *stmt_, const unsigned origin,
        const double pos, const double esize,
        const vector<pair<unsigned, double>> &freqs)
    {
        if (stmt_ == NULL)
            {
                throw std::runtime_error("NULL statement encountered");
            }
        for (auto &&fi : freqs)
            {
                int idx = 1;
                sqlite3_bind_int(
                    stmt_, idx++,
                    static_cast<int>(this->label)); // replicate ID
                sqlite3_bind_int(stmt_, idx++,
                                 static_cast<int>(fi.first)); // generation
                sqlite3_bind_int(stmt_, idx++, static_cast<int>(origin));
                sqlite3_bind_double(stmt_, idx++, pos);
                sqlite3_bind_double(stmt_, idx++, esize);
                sqlite3_bind_double(stmt_, idx++, fi.second); // freq
                int rc = sqlite3_step(stmt_);
                if (rc != SQLITE_DONE)
                    {
                        ostringstream message;
                        message << "Prepared statement could not be executed: "
                                << this_thread::get_id() << ' ' << __LINE__
                                << ' ' << idx << ' ' << rc << '\n';
                        throw std::runtime_error(message.str().c_str());
                    }
                sqlite3_reset(stmt_);
            }
        return static_cast<unsigned>(freqs.size());
    }

    void
    trajSQLonedb::call_operator_details(
        const fwdpy::selected_mut_tracker::final_t &data)
    {
        unsigned nrecords_passed = 0;

        for (auto &&outer : data)
            {
                if (this->tf->apply_origin_filter(outer.first))
                    {
                        for (auto &&inner : outer.second)
                            {
                                if (this->tf->apply_pos_esize_filter(
                                        inner.first)
                                    && this->tf->apply_freq_filter(
                                           inner.second))
                                    {
                                        nrecords_passed
                                            += apply_prepared_statement(
                                                memdb, memdb_stmt, outer.first,
                                                inner.first.first,
                                                inner.first.second,
                                                inner.second);
                                        if (nrecords_passed > threshold)
                                            {
                                                lock_guard<mutex> lock(
                                                    *dblock);
                                                copy_from_mem_db(
                                                    sql_callback_onedb);
                                                nrecords_passed = 0;
                                            }
                                    }
                            }
                    }
            }
        if (nrecords_passed)
            {
                lock_guard<mutex> lock(*dblock);
                copy_from_mem_db(sql_callback_onedb);
            }
    }

    // database columns are:
    // rep (int): only if onedb == true
    // generation (int)
    // origin (int)
    // pos (real)
    // esize (real)
    // freq(real)
    //
    // Indexes will be:
    // generation
    // rep: only if onedb == true
    //
    // The table name will be 'freqs'
    string
    traj2sql_details(const shared_ptr<mutex> &dblock_,
                     const fwdpy::selected_mut_tracker &data,
                     const fwdpy::trajFilter *tf, const string &dbname,
                     unsigned threshold, const unsigned label,
                     const bool onedb, const bool append)
    {
        try
            {
                if (onedb)
                    {
                        trajSQLonedb t(dblock_, tf, dbname, threshold, label,
                                       append);
                        t.prepare_statements();
                        t(data.final());
                    }
                else
                    {
                        ostringstream db;
                        db << dbname << '.' << label << ".db";
                        auto name = db.str();
                        trajSQL t(NULL, tf, name, threshold, append);
                        t.prepare_statements();
                        t(data.final());
                    }
            }
        catch (std::runtime_error &re)
            {
                return string(re.what());
            }
        return string();
    }
}

namespace fwdpy
{
    bool
    all_origins_pass(const unsigned)
    {
        return true;
    }
    bool
    all_pos_esize_pass(const std::pair<double, double> &)
    {
        return true;
    }
    bool
    all_freqs_pass(const std::vector<std::pair<KTfwd::uint_t, double>> &)
    {
        return true;
    }
    void
    traj2sql(const vector<unique_ptr<fwdpy::sampler_base>> &samplers,
             const shared_ptr<mutex> &dblock, const fwdpy::trajFilter *tf,
             const string &dbname, unsigned threshold, const unsigned label,
             const bool onedb, const bool append)
    {
        vector<future<string>> tasks;
        unsigned dummy = 0;
        for (auto &&i : samplers)
            {
                tasks.emplace_back(async(
                    launch::async, traj2sql_details, dblock,
                    *dynamic_cast<fwdpy::selected_mut_tracker *>(i.get()), tf,
                    dbname, threshold, label + dummy, onedb, append));
                dummy++;
            }

        ostringstream errors;
        for (auto &t : tasks)
            {
                auto m = t.get();
                if (!m.empty())
                    {
                        errors << m << '\n';
                    }
            }
        if (!errors.str().empty())
            {
                throw runtime_error(errors.str());
            }
    }
}
