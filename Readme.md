# Clone the repo and the submodules
```sh
git clone git@github.com:gbinux/ibdtools.git --recurse-submodules
cd ibdtools
```
# Dependencies
1. `htslib`
2. `fmt` 
3. `gtest`

They can be install using conda
```sh
conda env create -f env.yml
```
# Compile ibdtools

```sh
conda activate ibdtools
cd src
make ibdtools
```

# Simulated data and example

```sh
cd example

# simulating data
./simulate_data.py

# running ibdtools on the simulated input data
example/example.sh
```

# Documentation

1. List all available subcommand:
```sh
./ibdtools 
```

2. Show help message/documentation for a subcommand
```sh
./ibdtools [subcommand] -h
```
