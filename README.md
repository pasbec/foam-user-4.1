# README
------



## Installation
------

### Ubuntu 14.04

##### Pre-requirements
```bash
sudo apt-get install git-core build-essential binutils-dev cmake flex bison zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev libxt-dev rpm mercurial graphviz openmpi-bin python-all mayavi2
```

#### Basic

##### Download/Clone
```bash
cd "$HOME"
#
git clone --recursive 'ssh://git@git.becksteins.de:2222/pascal/foam-extend.git' 'foam'
```
##### Git user data
```bash
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
##### Git repository for dotfiles (optional)
```bash
cd "$HOME/foam"
#
git remote add 'dotfiles' 'ssh://git@git.becksteins.de:2222/pascal/dotfiles.git'
```

#### OpenFOAM Extend

##### Git repositories (optional)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
git remote add 'upstream' 'https://git.code.sf.net/p/openfoam-extend/foam-extend-3.1'
git remote add 'upstream_mirror' 'https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-foam-extend-3.1.git'
git remote add 'sf' 'ssh://p-be@git.code.sf.net/u/p-be/openfoam-extend'
#git remote set-url 'origin' 'https://git.becksteins.de/foam-extend/foam-extend-3.1'
```
##### Official documentation
```bash
cd "$HOME/foam/foam-extend-3.1"
#
less 'doc/buildInstructions/Ubuntu/Ubuntu_14.04'
```
##### Prepare ParaView compilation (qmake)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
export QT_BIN_DIR="$(dirname $(which qmake))"
echo "export QT_BIN_DIR=$(dirname "$(which qmake)")" >> 'etc/prefs.sh'
cd -
```
