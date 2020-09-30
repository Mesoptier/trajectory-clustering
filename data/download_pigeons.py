#!/usr/bin/env python3

import pathlib
import requests
import tarfile

def download():
    req = requests.get('https://royalsocietypublishing.org/action/'
        'downloadSupplement?doi=10.1098%2Frsbl.2013.0885&'
        'file=rsbl-2013-0885-file006.gz')
    with open('pigeon_data.tar', 'wb') as out:
        for chunk in req.iter_content(chunk_size=128):
            out.write(chunk)

def main():
    download()
    tarfile.open('pigeon_data.tar').extractall()
    pathlib.Path('pigeon_data.tar').unlink()

if __name__ == "__main__":
    main()
