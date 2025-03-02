#import library
import sqlite3

#create DatabaseManager class to create, load and run query against SQLite database
class DatabaseManager:
    #Create a constructor
    #df_file is the path to connect to SQL database
    def __init__(self, db_file):
        self.db_file = db_file
    
    #Function to create a database
    def createdb(self, sql_script):
        #Accepts DDL statements in the form of an sql script and creates a database
        #sql_script: the DDL statements stored in .sql file
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        error = False
        try:
            cur.executescript(sql_script)
        except Exception as e:
            print(f'Exception (DatabaseManager.py): {e}')
            error = True
            connection.rollback()
            raise SystemExit(1)
        if not error:
            connection.commit()
        cur.close()
        connection.close()
        return cur
    
    #Function to return query
    def query_db(self, sql, params):
        #Execeutes the SELECT statement for a given SQL database
        #sql: query for execution in the form of SELECT statement
        #param: Parameters for the query
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        error = False
        try:
            cur.execute(sql, params)
            rows = cur.fetchall()
        except Exception as e:
            print(f'Exception (DatabaseManager.py): {e}')
            error = True
            connection.rollback()
            raise SystemExit(1)
        if not error:
            connection.commit()
        cur.close()
        connection.close()
        return rows
        
    
    def insertdb(self, sql, params):
        #Executes the INSERT statements to load data into the SQL database 
        #sql: INSERT statement to load data into the database
        #param: Parameters for the INSERT statement.
        connection = sqlite3.connect(self.db_file)
        cur = connection.cursor()
        error = False
        try:
            cur.executemany(sql, params)
        except Exception as e:
            print(f'Exception (DatabaseManager.py): {e}')
            error = True
            connection.rollback()
            raise SystemExit(1)
        if not error:
            connection.commit()
        cur.close()
        connection.close()




