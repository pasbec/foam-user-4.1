# README
------

## Installation
------

### Pre-requirements

#### Ubuntu
##### Ubuntu 12.04 <a name="ubuntu1204"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev flex \
bison zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl
```
##### Ubuntu 13.10 <a name="ubuntu1310"></a> ([Ubuntu 12.04](#ubuntu1204))
##### Ubuntu 14.04 <a name="ubuntu1404"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
bison zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl
```
##### Ubuntu 15.10 <a name="ubuntu1510"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl
```
##### Ubuntu 16.04 <a name="ubuntu1604"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl

# For Paraview 5.4.1 with system qt
sudo apt-get install qt5-default libqt5x11extras5 libqt5x11extras5-dev
```

### Download/Clone root repository

```bash
cd "$HOME"
#
git clone --recursive 'ssh://git@git.becksteins.de:2222/pascal/foam-extend.git' 'foam'
```

### Git settings (optional)

#### User data
```bash
cd "$HOME/foam"
#
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
#### Add remote for dotfiles-repository
```bash
cd "$HOME/foam"
#
git remote add 'dotfiles' 'ssh://git@git.becksteins.de:2222/pascal/dotfiles.git'
```
#### Pull changes for bashrc.d subtree in user folder from dotfiles
```bash
cd "$HOME/foam"
#
git subtree pull -P 'pascal-4.0/etc/bashrc.d' 'dotfiles' 'subtree_dot/bashrc.d' --squash
```
#### Update submodule foam-extend
```bash
cd "$HOME/foam"
#
git submodule update --remote 'foam-extend-4.0'
```
#### User data for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-4.0"
#
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
#### Git repositories for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-4.0"
#
#git remote set-url 'origin' 'https://git.becksteins.de/foam-extend/foam-extend-4.0'
git remote add 'sf' 'ssh://p-be@git.code.sf.net/u/p-be/foam-extend-4.0'
git remote add 'upstream' 'https://git.code.sf.net/p/foam-extend/foam-extend-4.0'
git remote add 'upstream_mirror' 'https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0.git'
```


### Project part


#### Prepare preferences from template
```bash
cd "$HOME/foam/foam-extend-4.0"
#
# Template
cp 'etc/prefs.sh-EXAMPLE' 'etc/prefs.sh'
# Bison 2.7 (only Ubuntu 15.10+, maybe not needed any more)
#sed -i "s/#\(export WM_THIRD_PARTY_USE_BISON_27=1\)/\1/g" 'etc/prefs.sh'
#perl -pe "s{#(export WM_THIRD_PARTY_USE_BISON_27=1}{\1}" 'etc/prefs.sh'
# System qmake
#sed -i "/#export QT_BIN_DIR=\$QT_DIR\/bin/a export QT_BIN_DIR=$(dirname $(which qmake))" 'etc/prefs.sh'
perl -pe "s{#export QT_BIN_DIR=(\\$)QT_DIR/bin}{unset QT_THIRD_PARTY\nexport QT_BIN_DIR=$(dirname $(which qmake))}" 'etc/prefs.sh'
# Paraview 5.4.1
perl -pe "s{#(export WM_THIRD_PARTY_USE_CMAKE_332=1)}{unset WM_THIRD_PARTY_USE_CMAKE_322\n\1}" 'etc/prefs.sh'
perl -pe "s{#(export WM_THIRD_PARTY_USE_QT_580=1)}{unset WM_THIRD_PARTY_USE_QT_486\n\1}" 'etc/prefs.sh'
perl -pe "s{#(export WM_THIRD_PARTY_USE_PARAVIEW_541=1)}{unset WM_THIRD_PARTY_USE_PARAVIEW_440\n\1}" 'etc/prefs.sh'
cd -
```
#### Source environment
```bash
cd "$HOME/foam/foam-extend-4.0"
#
source 'etc/bashrc'
```
#### Install and configure ccache (optional)
```bash
sudo apt-get install ccache
export PATH="/usr/lib/ccache:${PATH}"
```
#### Compilation/Installation
```bash
cd "$HOME/foam/foam-extend-4.0"
#
./Allwmake.firstInstall | tee Allwmake.firstInstall.log
./Allwmake
```


### User part

#### Source environment
```bash
cd "$HOME/foam/pascal-4.0"
#
export FE40_FOAM_INST_DIR="$HOME/foam"
export FE40_FOAM_ETC="$FE40_FOAM_INST_DIR/foam-extend-4.0/etc"
export FE40_FOAM_USER_ETC="$HOME/foam/$USER-4.0/etc"
source "$FE40_FOAM_USER_ETC/bashrc"
fe40
```
#### Compilation/Installation
```bash
cd "$HOME/foam/pascal-4.0"
#
./Allwmake.firstInstall | tee Allwmake.firstInstall.log
./Allwmake
```
