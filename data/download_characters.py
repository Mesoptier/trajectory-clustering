#!/usr/bin/env python3

import pathlib
import requests
import scipy.io

def download():
    req = requests.get('https://archive.ics.uci.edu/ml/'
        'machine-learning-databases/character-trajectories/'
        'mixoutALL_shifted.mat')
    with open('mixoutALL_shifted.mat', 'wb') as out:
        for chunk in req.iter_content(chunk_size=128):
            out.write(chunk)

def main():
    pathlib.Path('characters/data').mkdir(parents=True, exist_ok=True)
    download()
    mat = scipy.io.loadmat('mixoutALL_shifted.mat')
    n = 2858
    with open('characters/data/dataset.txt', 'w') as dataset:
        for i in range(n):
            index = mat['consts']['charlabels'][0][0][0][i] - 1
            label = mat['consts']['key'][0][0][0][index][0]
            xs = mat['mixout'][0][i][0]
            ys = mat['mixout'][0][i][1]
            name = label + '{:04d}'.format(i + 1) + '.txt'
            print(name, file=dataset)
            with open('characters/data/' + name, 'w') as out:
                for (x, y) in zip(xs, ys):
                    print(x, y, file=out)
    pathlib.Path('mixoutALL_shifted.mat').unlink()

if __name__ == "__main__":
    main()
