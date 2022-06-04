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
