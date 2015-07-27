# README
------

## Installation
------

### Ubuntu 14.04

#### Basic

##### Pre-requirements
```bash
sudo apt-get install git-core build-essential binutils-dev cmake flex bison \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev libxt-dev \
rpm mercurial graphviz openmpi-bin python-all mayavi2
```
##### Download/Clone
```bash
cd "$HOME"
#
git clone --recursive 'ssh://git@git.becksteins.de:2222/pascal/foam-extend.git' 'foam'
```

#### Git settings (optional)

##### User data
```bash
cd "$HOME/foam"
#
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
##### Add remote for dotfiles-repository
```bash
cd "$HOME/foam"
#
git remote add 'dotfiles' 'ssh://git@git.becksteins.de:2222/pascal/dotfiles.git'
```
##### Pull changes for bashrc.d subtree in user folder from dotfiles
```bash
cd "$HOME/foam"
#
git subtree pull -P 'pascal-3.1/etc/bashrc.d' 'dotfiles' 'subtree_dot-bashrc.d' --squash
```
##### Update submodule foam-extend
```bash
cd "$HOME/foam"
#
git submodule update --remote 'foam-extend-3.1'
```
##### User data for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
##### Git repositories for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
git remote add 'upstream' 'https://git.code.sf.net/p/openfoam-extend/foam-extend-3.1'
git remote add 'upstream_mirror' 'https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-foam-extend-3.1.git'
git remote add 'sf' 'ssh://p-be@git.code.sf.net/u/p-be/openfoam-extend'
#git remote set-url 'origin' 'https://git.becksteins.de/foam-extend/foam-extend-3.1'
```


#### Project: OpenFOAM Extend

##### Official documentation
```bash
cd "$HOME/foam/foam-extend-3.1"
#
less 'doc/buildInstructions/Ubuntu/Ubuntu_14.04'
```
##### Source environment
```bash
cd "$HOME/foam/foam-extend-3.1"
#
source 'etc/bashrc'
```
##### Prepare ParaView compilation (qmake)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
export QT_BIN_DIR="$(dirname $(which qmake))"
echo "export QT_BIN_DIR=$(dirname "$(which qmake)")" >> 'etc/prefs.sh'
cd -
```
##### Installation
```bash
cd "$HOME/foam/foam-extend-3.1"
#
./Allwmake.firstInstall
```


#### User

##### Source environment
```bash
cd "$HOME/foam/pascal-3.1"
#
source 'etc/bashrc'
fe31
```