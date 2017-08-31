import os, json, itertools

FILEPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
FILENAME = os.path.join(FILEPATH, 'go.json')

class Searcher(object):

    def __init__(self):
        with open(FILENAME, 'r') as infile:
            self.go = json.load(infile)

    def traverse(self, term, inclusive=True):
        if inclusive:
            relatives = [term]
        else:
            relatives = []

        new_relatives = self.go[term]['children']
        print(new_relatives)

        while True:
            update = list(set(relatives).union(new_relatives))
            if len(relatives) == len(update):
                break
            else:
                relatives = update
            new_relatives = itertools.chain.from_iterable([self.go[t]['children'] for t in new_relatives])

        return relatives

    def select_namespace(self, namespace):
        return [term for term in self.go if self.go[term]['namespace'] == namespace]

