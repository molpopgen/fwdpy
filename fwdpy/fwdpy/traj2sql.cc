#include <thread>
#include <mutex>
#include <functional>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <cstdio>
#include <iostream>
#include <sqlite3.h>
#include "sampler_selected_mut_tracker.hpp"

using namespace std;

namespace
{
    // when onedb == true, we
    // need to synchronise write
    // access to our db
    mutex dblock;

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
        rc = sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL,
                          &error_message);
        return rc;
    }

    int
    execute_sql_statement(sqlite3 *db, const string &sql, char *error_message)
    {
        return sqlite3_exec(db, sql.c_str(), NULL, NULL, &error_message);
    }

    // struct user_data_wrapper
    //{
    //    sqlite3 *recipient;
    //    sqlite3_stmt *recipient_stmt;
    //    explicit user_data_wrapper(sqlite3 *db_, sqlite3_stmt *stmt_)
    //        : recipient(db_), recipient_stmt(stmt_)
    //    {
    //    }
    //};

    static int
    sql_callback(void *data, int argc, char **argv, char **azColName)
    {
		//std::cout << "about 2 copy\n";
		//for(int i=0;i<argc;++i)
		//{
		//	cout << argv[i] << ' ';
		//}
		//cout << '\n';
        auto dp = reinterpret_cast<sqlite3_stmt *>(data);
        int idx = 1;
        sqlite3_bind_int(dp, idx, atoi(argv[idx-1])); // generation
        idx++;
        sqlite3_bind_int(dp, idx, atoi(argv[idx-1]));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx-1], NULL));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx-1], NULL));
        idx++;
        sqlite3_bind_double(dp, idx, strtod(argv[idx-1], NULL)); // freq
		int rc = sqlite3_step(dp);
		if(rc!=SQLITE_DONE)
		{
			std::cout << "oh crap\n";
		}
        sqlite3_reset(dp);
        return 0;
    };

    class trajSQL
    {
      private:
        sqlite3 *db, *memdb;
        sqlite3_stmt *stmt, *memdb_stmt;
        char *error_message;
        unsigned threshold;
      public:
        using callback_fxn_t = int (*)(void *, int, char **, char **);
        explicit trajSQL(sqlite3 *db_, const string &dbname_,
                         const unsigned threshold_, const bool append_);
        virtual ~trajSQL() noexcept;
        void handle_return_values(int rc, int line_num);
        virtual void create_table(sqlite3 *db_);
        virtual void create_index(sqlite3 *db_);
        virtual void prepare_statements();
        virtual void copy_from_mem_db(callback_fxn_t);
        unsigned
        apply_prepared_statement(sqlite3 *db_, sqlite3_stmt *stmt_,
                                 const unsigned origin, const unsigned pos,
                                 const unsigned esize,
                                 const vector<pair<unsigned, double>> &freqs);
        virtual void
        operator()(const fwdpy::selected_mut_tracker::final_t &data);
    };

    trajSQL::trajSQL(sqlite3 *db_, const string &dbname_,
                     const unsigned threshold_, const bool append_)
        : db(nullptr), memdb(nullptr), stmt(nullptr), memdb_stmt(nullptr),
          error_message(nullptr), threshold(threshold_)
    {
        //cout << "in constructor" << endl;
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
                //cout << "opened db " << dbname_ << '\n';
            }

        int rc = sqlite3_open(":memory:", &memdb);
        handle_return_values(rc, __LINE__);
        create_table(memdb);
        //cout << "opened :memory:" << endl;
        rc = apply_sql_pragma(memdb, error_message);
        handle_return_values(rc, __LINE__);
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
        //cout << "handling a return value : " << this_thread::get_id() << ' '
        //     << (rc == SQLITE_OK) << ' ' << line_num << ' '
        //     << (error_message == NULL) << ' ' << (error_message == nullptr)
        //     << '\n';
        if (rc != SQLITE_OK)
            {
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
		string prepped_statement="insert into freqs values (?1,?2,?3,?4,?5);";
        sqlite3_prepare_v2(db, prepped_statement.c_str(),
                           prepped_statement.size(), &stmt, NULL);
		if(stmt==NULL)
		{
			throw runtime_error("could not prepare statement for on-disk database.");
		}
        sqlite3_prepare_v2(memdb, prepped_statement.c_str(),
                           prepped_statement.size(), &memdb_stmt, NULL);
		if(memdb_stmt==NULL)
		{
			throw runtime_error("could not prepare statement for in-memory database.");
		}
    }


    unsigned
    trajSQL::apply_prepared_statement(
        sqlite3 *db_, sqlite3_stmt *stmt_, const unsigned origin,
        const unsigned pos, const unsigned esize,
        const vector<pair<unsigned, double>> &freqs)
    {
		if(stmt_ == NULL)
		{
			throw std::runtime_error("NULL statement encountered");
		}
        for (auto &&fi : freqs)
            {
                int idx = 1;
                //cout << "insertions: " << this_thread::get_id() << ' '
                //     << fi.first << ' ' << origin << ' ' << pos << ' ' << esize
                //     << ' ' << fi.second << endl;
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
							<<this_thread::get_id() << ' ' << __LINE__ << ' ' << idx << ' ' << rc << '\n';
                        throw std::runtime_error(message.str().c_str());
                    }
                sqlite3_reset(stmt_);
            }
        return static_cast<unsigned>(freqs.size());
    }

    void trajSQL::copy_from_mem_db(callback_fxn_t)
    {
        string sql = "select * from freqs";
		int rc=sqlite3_exec(db,"BEGIN TRANSACTION;",NULL,NULL,&error_message);
		handle_return_values(rc,__LINE__);
        rc = sqlite3_exec(memdb, sql.c_str(), sql_callback, (void *)stmt,
                              &error_message);
        handle_return_values(rc, __LINE__);
		rc=sqlite3_exec(db,"COMMIT TRANSACTION;",NULL,NULL,&error_message);
		handle_return_values(rc,__LINE__);
        sql = "delete from freqs;";
        rc = sqlite3_exec(memdb, sql.c_str(), NULL, NULL, &error_message);
        handle_return_values(rc, __LINE__);
    }

    void
    trajSQL::operator()(const fwdpy::selected_mut_tracker::final_t &data)
    {
        unsigned nrecords_passed = 0;
        for (auto &&outer : data)
            {
                for (auto &&inner : outer.second)
                    {
                        nrecords_passed += apply_prepared_statement(
                            memdb, memdb_stmt, outer.first, inner.first.first,
                            inner.first.second, inner.second);
                        if (nrecords_passed >= threshold)
                            {
                                copy_from_mem_db(sql_callback);
                                threshold = 0;
                            }
                    }
            }
        if (nrecords_passed)
            {
                copy_from_mem_db(sql_callback);
            }
    }
    //
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
    void
    traj2sql_details(const fwdpy::selected_mut_tracker &data,
                     fwdpy::origin_filter_fxn origin_filter,
                     fwdpy::pos_esize_filter_fxn pos_esize_filter,
                     fwdpy::freq_filter_fxn freq_filter, const string &dbname,
                     unsigned threshold, const unsigned label,
                     const bool onedb, const bool append)
    {
        // trajSQL(sqlite3 *db_, const string &dbname_,
        //                 const unsigned threshold_, const bool append_);
        if (onedb)
            {
            }
        else
            {
                try
                    {
                        ostringstream db;
                        db << dbname << '.' << label << ".db";
                        auto name = db.str();
                        //cout << name << '\n';
                        trajSQL t(NULL, name, threshold, append);
						t.prepare_statements();
                        t(data.final());
                    }
                catch (std::runtime_error &re)
                    {
                        throw re;
                    }
            }
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
             fwdpy::origin_filter_fxn origin_filter,
             fwdpy::pos_esize_filter_fxn pos_esize_filter,
             fwdpy::freq_filter_fxn freq_filter, const string &dbname,
             unsigned threshold, const unsigned label, const bool onedb,
             const bool append)
    {
        vector<thread> threads;

        unsigned dummy = 0;
        for (auto &&i : samplers)
            {
                threads.emplace_back(
                    traj2sql_details,
                    *dynamic_cast<fwdpy::selected_mut_tracker *>(i.get()),
                    origin_filter, pos_esize_filter, freq_filter, dbname,
                    threshold, label + dummy, onedb, append);
                dummy++;
            }

        for (auto &t : threads)
            {
                t.join();
            }
    }
}
