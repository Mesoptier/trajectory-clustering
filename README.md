This is the repository for the code that supplements the paper
“(k, ℓ)-Medians Clustering of Trajectories Using Continuous Dynamic Time Warping”
by Milutin Brankovic, Kevin Buchin, Koen Klaren, André Nusser, Aleksandr Popov,
and Sampson Wong, published in _28th International Conferenceon Advances in
Geographic Information Systems (SIGSPATIAL ’20)._

# Prerequisites
Make sure you have a modern C++ compiler that supports C++ 17.
GCC, clang, or MSVC from the last few years should be fine.
Make sure you also have Python 3 with scipy and matplotlib installed to get the
data and do the plotting.
Other tools needed: git, cmake, make.

# Setting up
To perform the experiments yourself, do the following steps.
1. Clone the latest state of the master branch:
```
git clone --depth 1 --recurse-submodules --shallow-submodules https://github.com/Mesoptier/trajectory-clustering.git clustering && cd clustering
```
2. Get the character data.
```
cd data && ./download_characters.py
```
3. Get the pigeon data.
```
./download_pigeons.py
./coord_conversion.R
```
4. Get the stork data.
```
./download_movebank.py
```
5. Build the project.
```
mkdir ../build && cd ../build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j6
cd ..
```

# The experiments
To run the experiments, simply execute `./build/TrajectoryClustering_run`.
There is some extra explanation printed out to stdout when running the program.
