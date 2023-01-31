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
sudo apt-get install freecad
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
sudo apt-get install freecad
```
##### Ubuntu 15.10 <a name="ubuntu1510"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl
sudo apt-get install freecad
```
##### Ubuntu 16.04 <a name="ubuntu1604"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl
sudo apt-get install freecad

# For Paraview 5.4.1 with system qt
sudo apt-get install qt5-default libqt5x11extras5 libqt5x11extras5-dev
```
##### WSL Ubuntu 16.04 <a name="ubuntu1604"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin libopenmpi-dev python-all mayavi2
sudo apt-get install perl
sudo apt-get install freecad

# To avoid problems
sudo apt-get install mercurial subversion python-scipy lua5.2

# For Paraview 5.4.1 with system qt
sudo apt-get install qt5-default libqt5x11extras5 libqt5x11extras5-dev

# Modify foam-extend-4.1/ThirdParty/AllMake.stage5 - after line 86 add:
                hg revert --all --rev 3360
```

### Download/Clone root repository

```bash
cd "$HOME"
#
git clone --recursive 'ssh://git@git.becksteins.de:2222/pascal/foam-extend.git' 'foam'
```

### Project part


#### Prepare preferences from template
```bash
cd "$HOME/foam/foam-extend-4.1"
#
# Template
cp 'etc/prefs.sh-EXAMPLE' 'etc/prefs.sh'
# System qt
perl -i -pe "s{^(export QT_THIRD_PARTY=1)}{#\1}" 'etc/prefs.sh'
perl -i -pe "s{(#export QT_BIN_DIR=(\\$)QT_DIR/bin)}{\1\nunset QT_THIRD_PARTY\nexport QT_BIN_DIR=$(dirname $(which qmake))}" 'etc/prefs.sh'
# Thirparty qt
perl -i -pe "s{#(export QT_THIRD_PARTY=1)}{\1}" 'etc/prefs.sh'
# Paraview 5.4.1
perl -i -pe "s{^(export WM_THIRD_PARTY_USE_CMAKE_322=1)}{#\1}" 'etc/prefs.sh'
perl -i -pe "s{#(export WM_THIRD_PARTY_USE_CMAKE_332=1)}{unset WM_THIRD_PARTY_USE_CMAKE_322\n\1}" 'etc/prefs.sh'
perl -i -pe "s{^(export WM_THIRD_PARTY_USE_QT_486=1)}{#\1}" 'etc/prefs.sh'
perl -i -pe "s{#(export WM_THIRD_PARTY_USE_QT_580=1)}{unset WM_THIRD_PARTY_USE_QT_486\n\1}" 'etc/prefs.sh'
perl -i -pe "s{^(export WM_THIRD_PARTY_USE_PARAVIEW_440=1)}{#\1}" 'etc/prefs.sh'
perl -i -pe "s{#(export WM_THIRD_PARTY_USE_PARAVIEW_541=1)}{unset WM_THIRD_PARTY_USE_PARAVIEW_440\n\1}" 'etc/prefs.sh'
cd -
```
#### Source environment
```bash
cd "$HOME/foam/foam-extend-4.1"
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
cd "$HOME/foam/foam-extend-4.1"
#
./Allwmake.firstInstall | tee 'Allwmake.firstInstall.log'
./Allwmake
```


### User part

#### Source environment
```bash
cd "$HOME/foam"
ln -s 'user-4.1' "$USER-4.1"
cd "$USER-4.1"
#
export FE41_FOAM_INST_DIR="$HOME/foam"
export FE41_FOAM_ETC="$FE41_FOAM_INST_DIR/foam-extend-4.1/etc"
export FE41_FOAM_USER_ETC="$HOME/foam/$USER-4.1/etc"
source "$FE41_FOAM_USER_ETC/bashrc"
fe41
```
#### Compilation/Installation
```bash
cd "$HOME/foam/$USER-4.1"
#
./Allwmake.firstInstall | tee 'Allwmake.firstInstall.log'
./Allwmake
```
