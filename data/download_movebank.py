#!/usr/bin/env python3

import pathlib
import requests

def download():
    req = requests.get('https://www.datarepository.movebank.org/'
        'bitstream/handle/10255/move.792/Eastern%20flyway%20spring%20migration'
        '%20of%20adult%20white%20storks%20%28data%20from%20Rotics%20et%20al.'
        '%202018%29-gps.csv')
    with open('storks.csv', 'wb') as out:
        for chunk in req.iter_content(chunk_size=128):
            out.write(chunk)

def main():
    pathlib.Path('movebank').mkdir(exist_ok=True)
    download()

    trajectories = {}
    with open('storks.csv', 'r') as infile:
        next(infile)
        for line in infile:
            if not line or line.isspace():
                continue
            point = line.split(',')
            year = point[2][:4]
            ident = point[13][1:-1]
            if ident + '_' + year not in trajectories:
                trajectories[ident + '_' + year] = []
            trajectories[ident + '_' + year].append(point[3:5])

    with open('movebank/dataset.txt', 'w') as ds:
        for key in trajectories:
            print(key + '.csv', file=ds)
            with open('movebank/' + key + '.csv', 'w') as tout:
                print('x, y', file=tout)
                for x, y in trajectories[key]:
                    print(x, y, sep=', ', file=tout)

    pathlib.Path('storks.csv').unlink()

if __name__ == '__main__':
    main()
