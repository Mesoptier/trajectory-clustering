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
    base = pathlib.Path('Data_for_Mann_et_al_RSBL')
    (base / 'Horspath' / 'H71' / 'H71_18_16_08_07.txt').unlink()
    (base / '40Kmap.jpg').unlink()
    for f in base.rglob('*.DS_Store'):
        f.unlink()

    subdirs = ['Bladon & Church route recapping/bladon heath',
        'Bladon & Church route recapping/church hanborough',
        'Horspath', 'Weston']
    for sd in subdirs:
        p = base / sd
        with open(p / 'dataset.txt', 'w') as ds:
            for x in sorted(p.rglob('*.txt')):
                if x.name != 'dataset.txt' and 'Untrimmed' not in x.parts and \
                        'untrimmed' not in x.parts:
                    print(x.relative_to(p), file=ds)
        pigeon_dirs = [d for d in p.iterdir() if d.is_dir()]
        for psd in pigeon_dirs:
            with open(psd / 'dataset.txt', 'w') as ds:
                for x in sorted((psd).glob('*.txt')):
                    if x.name != 'dataset.txt':
                        print(x.name, file=ds)

if __name__ == '__main__':
    main()
