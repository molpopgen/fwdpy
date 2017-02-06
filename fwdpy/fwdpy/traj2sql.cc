#include <thread>
#include <mutex>
#include <functional>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <sqlite3.h>
#include "sampler_selected_mut_tracker.hpp"

using namespace std;

namespace
{
    // when onedb == true, we
    // need to synchronise write
    // access to our db
    mutex dblock;

    void
    apply_sql_pragma(sqlite3 *db)
    {
        char *error_message;
        // from http://blog.quibb.org/2010/08/fast-bulk-inserts-into-sqlite/
        int rc = sqlite3_exec(db, "PRAGMA synchronous=OFF", NULL, NULL,
                              &error_message);

        if (rc != SQLITE_OK)
            {
                string error(error_message);
                sqlite3_free(error_message);
                throw std::runtime_error(error.c_str());
            }
        rc = sqlite3_exec(db, "PRAGMA count_changes=OFF", NULL, NULL,
                          &error_message);
        if (rc != SQLITE_OK)
            {
                string error(error_message);
                sqlite3_free(error_message);
                throw std::runtime_error(error.c_str());
            }
        rc = sqlite3_exec(db, "PRAGMA journal_mode=MEMORY", NULL, NULL,
                          &error_message);
        if (rc != SQLITE_OK)
            {
                string error(error_message);
                sqlite3_free(error_message);
                throw std::runtime_error(error.c_str());
            }
        rc = sqlite3_exec(db, "PRAGMA temp_store=MEMORY", NULL, NULL,
                          &error_message);
        if (rc != SQLITE_OK)
            {
                string error(error_message);
                sqlite3_free(error_message);
                throw std::runtime_error(error.c_str());
            }
    }

	void create_index(sqlite3 * db,const bool onedb)
	{
		string sql;
		if(onedb)
		{
			sql = "create index rep_gen on freqs(rep,generation);";
		}
		else
		{
			sql = "create index gen on freqs(generation);";
		}
		char * error_message;
		int rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&error_message);
		if(rc!=SQLITE_OK)
		{
			string error(error_message);
			sqlite3_free(error_message);
			throw std::runtime_error(error.c_str());
		}
	}

    sqlite3 *
    open_database(const string &dbname, const unsigned label, const bool onedb)
    {
        sqlite3 *db;

        if (onedb == false)
            {
                ostringstream n;
                n << dbname << '.' << label << ".db";
                int i = sqlite3_open(n.str().c_str(), &db);
                if (i)
                    {
                        throw std::runtime_error("could not open database "
                                                 + n.str());
                    }
				apply_sql_pragma(db);
            }
        else
            {
                int i = sqlite3_open(dbname.c_str(), &db);
                if (i)
                    {
                        throw std::runtime_error("could not open database "
                                                 + dbname);
                    }
				lock_guard<mutex> lock(dblock);
				apply_sql_pragma(db);
            }
        return db;
    }

    int
    sql_callback(void *NotUsed, int argc, char **argv, char **azColName)
    {
        return 0;
    };

    void
    execute_sql_statement(sqlite3 *db, const string &sql, const bool onedb)
    {
        char *error_message;
        if (onedb)
            {
                lock_guard<mutex> lock(dblock);
                int rc = sqlite3_exec(db, sql.c_str(), sql_callback, 0,
                                      &error_message);
                if (rc != SQLITE_OK)
                    {
                        string error(error_message);
                        sqlite3_free(error_message);
                        throw std::runtime_error(error.c_str());
                    }
            }
        else
            {
                int rc = sqlite3_exec(db, sql.c_str(), sql_callback, 0,
                                      &error_message);
                if (rc != SQLITE_OK)
                    {
                        string error(error_message);
                        sqlite3_free(error_message);
                        throw std::runtime_error(error.c_str());
                    }
            }
    }

    void
    create_table(sqlite3 *db, const bool onedb)
    {
        string sql("CREATE TABLE IF NOT EXISTS freqs(");
        if (onedb)
            {
                sql += "rep int NOT NULL,";
            }
        sql += "generation int NOT NULL, origin int NOT NULL, pos real NOT ";
        sql += "NULL, esize real NOT NULL, freq real NOT NULL);";
        execute_sql_statement(db, sql, onedb);
    }

    void
    initialize_insert_statement(ostringstream &sql, const bool onedb)
    {
        sql << "INSERT INTO FREQS (";
        if (onedb)
            {
                sql << "rep,";
            }
        sql << "generation,origin,pos,esize,freq) VALUES(";
    }

    unsigned
    update_insert_statement(ostringstream &sql, const bool onedb,
                            const unsigned label,
                            const KTfwd::uint_t origin_time, const double pos,
                            const double esize,
                            const vector<pair<KTfwd::uint_t, double>> &freqs)
    {
        for (auto &&fi : freqs)
            {
                initialize_insert_statement(sql, onedb);
                if (onedb)
                    {
                        sql << label << ',';
                    }
                sql << fi.first << ',' << origin_time << ',' << pos << ','
                    << esize << ',' << fi.second << ");\n";
            }
        return static_cast<unsigned>(freqs.size());
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
        sqlite3 *db = open_database(dbname, label, onedb);
        create_table(db, onedb);

        unsigned nrecords_passed = 0;
        ostringstream sql;
        for (auto &&outer : data)
            {
                if (origin_filter(outer.first))
                    {
                        for (auto inner : outer.second)
                            {
                                if (pos_esize_filter(inner.first))
                                    {
                                        if (freq_filter(inner.second))
                                            {
                                                nrecords_passed
                                                    += update_insert_statement(
                                                        sql, onedb, label,
                                                        outer.first,
                                                        inner.first.first,
                                                        inner.first.second,
                                                        inner.second);
                                                if (nrecords_passed
                                                    > threshold)
                                                    {
                                                        execute_sql_statement(
                                                            db, sql.str(),
                                                            onedb);
                                                        // reset buffer
                                                        sql.str(string());
                                                        nrecords_passed = 0;
                                                    }
                                            }
                                    }
                            }
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
