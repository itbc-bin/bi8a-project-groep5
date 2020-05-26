import mysql.connector


def connection_database():
    """Make connection to the database
    :return: The connection
    """
    connection = mysql.connector.connect(
        host='hannl-hlo-bioinformatica-mysqlsrv.mysql.database.azure.com',
        user='owe7_pg5@hannl-hlo-bioinformatica-mysqlsrv',
        database='owe7_pg5',
        password='blaat1234')

    return connection


def mand():
    connction = connection_database()
    cursor = connction.cursor()
    cursor.execute("select combination, amount from results join "
                   "algorithm_results ar on results.id = ar.results_id "
                   "where url_id = '{}'".format('mandje'))
    for result in cursor.fetchall():
        print(result[0], result[1])


if __name__ == '__main__':
    mand()
