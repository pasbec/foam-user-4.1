# README

## Git settings

### User data
```bash
git config user.name "Pascal Beckstein"
git config user.email "pascal@becksteins.de"
```

### Remotes

##### Base repository
```bash
git remote set-url origin ssh://git@git.becksteins.de:2222/pascal/foam-extend.git
```
##### Dependent repositories
```bash
git remote add dotfiles ssh://git@git.becksteins.de:2222/pascal/dotfiles.git
```
##### Submodule foam-extend in "foam-extend-3.1"
```bash
cd foam-extend-3.1
git remote add upstream https://git.code.sf.net/p/openfoam-extend/foam-extend-3.1
git remote add upstream_mirror https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-foam-extend-3.1.git
git remote add sf ssh://p-be@git.code.sf.net/u/p-be/openfoam-extend
git remote set-url origin https://git.becksteins.de/foam-extend/foam-extend-3.1
cd -
```

## Installation

### Ubuntu 14.04

##### Pre-requirements
```bash
sudo apt-get install git-core build-essential binutils-dev cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev libxt-dev rpm mercurial graphviz openmpi-bin python-all mayavi2
```



Default install of foem-extend-3.1 in $HOME/foam:
```bash
cd $HOME
mkdir foam
cd foam
git clone git://git.code.sf.net/p/openfoam-extend/foam-extend-3.1 foam-extend-3.1
```

##### ParaView (qmake)
```bash
export QT_BIN_DIR=$(dirname $(which qmake))
echo "export QT_BIN_DIR=$(dirname $(which qmake))" >> etc/prefs.sh
```
