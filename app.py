from Bio import Entrez
from flask import Flask, render_template, request

app = Flask(__name__)


# TODO
# return dictionary met search_term+gensymbool
# return lijst met alle sites ( searchterm+gensymbool)
# return
# functie om verbinding te maken met de database ( local database )

@app.route('/')
def hello_world():
    return render_template('index.html')


if __name__ == '__main__':
    app.run()
