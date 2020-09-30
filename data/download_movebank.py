#!/usr/bin/env python3

import pathlib
import requests

def download():
    req = requests.get('https://www.datarepository.movebank.org/'
        'bitstream/handle/10255/move.792/Eastern%20flyway%20spring%20migration'
        '%20of%20adult%20white%20storks%20%28data%20from%20Rotics%20et%20al.'
        '%202018%29-gps.csv')
    with open('storks.csv', 'w') as out:
        for chunk in req.iter_content(chunk_size=128):
            out.write(chunk)

def main():
    pathlib.Path('movebank/projection').mkdir(parents=True, exist_ok=True)
    download()
    # TODO
    pathlib.Path('storks.csv').unlink()

if __name__ == "__main__":
    main()
