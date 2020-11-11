#!/usr/bin/env python3

import pathlib
import subprocess

sites = [['a55', 'brc', 'c17', 'c35', 'p29', 'p39', 'p94'],
    ['a94', 'c22', 'c70', 'k77', 'l29', 'liv', 'r47', 's93'],
    ['H22', 'H27', 'H30', 'H35', 'H38', 'H41', 'H42', 'H71'],
    ['H23', 'H31', 'H32', 'H34', 'H36', 'H50', 'H58', 'H62']]
paths = ['Bladon & Church route recapping/bladon heath',
    'Bladon & Church route recapping/church hanborough',
    'Horspath', 'Weston']
names = ['bladon.gpx', 'church.gpx', 'horspath.gpx', 'weston.gpx']
base = pathlib.Path('convergence')
temp = base / 'tmp'

def rm_tree_dir(d) :
    for sub in d.iterdir():
        if sub.is_dir():
            rm_tree_dir(sub)
        else:
            sub.unlink()
    d.rmdir()

def prepare():
    temp.mkdir(exist_ok=True)
    for (i, p) in enumerate(paths):
        cur_tmp = temp / paths[i]
        cur_tmp.mkdir(parents=True, exist_ok=True)
        for s in sites[i]:
            all_pigeons_d = sorted((base / p).rglob(s + '_??_5'))
            if not all_pigeons_d:
                continue
            center = []
            with open(all_pigeons_d[-1], 'r') as file:
                file.readline()
                number_of_curves = int(file.readline())
                for j in range(number_of_curves):
                    file.readline()
                center = list(map(float, file.readline().rstrip().split(' ')))
            with open(cur_tmp / s, 'w') as tmpfile:
                print('x, y', file=tmpfile)
                for x, y in zip(center[0::2], center[1::2]):
                    print(x, y, sep=', ', file=tmpfile)

def write_gpx():
    out = base / 'gpx'
    out.mkdir(exist_ok=True)
    for (i, p) in enumerate(paths):
        cur_tmp = temp / paths[i]
        with open(out / names[i], 'w') as gpx:
            print('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                file=gpx)
            print('<gpx version="1.1">', file=gpx)
            for s in sites[i]:
                print('<trk>\n<trkseg>', file=gpx)
                with open(cur_tmp / s, 'r') as longlat:
                    longlat.readline()
                    for line in longlat:
                        lng, lat = line.split(',')
                        print('<trkpt lat="{0}" lon="{1}"></trkpt>'.format(
                            lat.strip(), lng.strip()), file=gpx)
                print('</trkseg>\n</trk>', file=gpx)
            print('</gpx>', file=gpx)

def main():
    prepare()
    subprocess.run("./coord_unconversion.R")
    write_gpx()
    rm_tree_dir(temp)

if __name__ == '__main__':
    main()
