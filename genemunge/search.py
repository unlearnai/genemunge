import os, json

FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
FILENAME = os.path.join(FILEPATH, 'go.json')


class Searcher(object):

    def __init__(self):
        with open(FILENAME, 'r') as infile:
            self.go = json.load(infile)
