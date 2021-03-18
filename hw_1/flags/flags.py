import argparse

class Args:
    pass

class ArgsParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='silly td-idf searcher')

        self.parser.add_argument(
                '-q', '--query',
                required = True,
                type = str,
                help = 'query string',
        )

        self.parser.add_argument(
                '-n', '--number',
                default=1,
                type=int,
                help='number of sentences relevant to the query'
        ) 
        
        self.parser.add_argument(
                '-d', '--download',
                default = False,
                type=bool,
                help='download text from links'
        )

        self.parser.add_argument(
                '-l', '--links',
                required = True,
                type = str,
                help = 'path to links file',
        )

        self.parser.add_argument(
                '-m', '--method',
                default = 1,
                type = int,
                help = 'method for tf',
        )
    
    def Parse(self, args):
        a = Args()
        return self.parser.parse_args(args=args, namespace=a)

