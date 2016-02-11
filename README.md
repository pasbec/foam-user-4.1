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
sudo apt-get install openmpi-bin python-all mayavi2
```
##### Ubuntu 13.10 <a name="ubuntu1310"></a> ([Ubuntu 12.04](#ubuntu1204))
##### Ubuntu 14.04 <a name="ubuntu1404"></a>
```bash
sudo apt-get update
sudo apt-get install git-core build-essential binutils-dev cmake flex \
bison zlib1g-dev qt4-dev-tools libqt4-dev libncurses5-dev libiberty-dev \
libxt-dev rpm mercurial graphviz
sudo apt-get install openmpi-bin python-all mayavi2
```
##### Ubuntu 15.10 <a name="ubuntu1510"></a> ([Ubuntu 14.04](#ubuntu1404))

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
git subtree pull -P 'pascal-3.2/etc/bashrc.d' 'dotfiles' 'subtree_dot/bashrc.d' --squash
```
#### Update submodule foam-extend
```bash
cd "$HOME/foam"
#
git submodule update --remote 'foam-extend-3.2'
```
#### User data for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-3.2"
#
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
#### Git repositories for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-3.2"
#
git remote add 'upstream' 'git://git.code.sf.net/p/foam-extend/foam-extend-3.2'
git remote add 'upstream_mirror' 'https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-3.2.git'
git remote add 'sf' 'https://p-be@git.code.sf.net/u/p-be/foam-extend-3.2'
#git remote set-url 'origin' 'https://git.becksteins.de/foam-extend/foam-extend-3.2'
```


### Project part


#### Prepare preferences from template
```bash
cd "$HOME/foam/foam-extend-3.2"
#
# Template
cp etc/prefs.sh-EXAMPLE etc/prefs.sh
# Bison 2.7
sed -i "s/#\(export WM_THIRD_PARTY_USE_BISON_27=1\)/\1/g" 'etc/prefs.sh'
# System qmake
sed -i "/#export QT_BIN_DIR=\$QT_DIR\/bin/a export QT_BIN_DIR=$(dirname $(which qmake))" 'etc/prefs.sh'
cd -
```
#### Source environment
```bash
cd "$HOME/foam/foam-extend-3.2"
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
cd "$HOME/foam/foam-extend-3.2"
#
./Allwmake.firstInstall
```


### User part

#### Source environment
```bash
cd "$HOME/foam/pascal-3.2"
#
export FE32_FOAM_INST_DIR="$HOME/foam"
export FE32_FOAM_ETC="$FE31_FOAM_INST_DIR/foam-extend-3.2/etc"
export FE32_FOAM_USER_ETC="$HOME/foam/$USER-3.2/etc"
source "$FE32_FOAM_USER_ETC/bashrc"
fe32
```
#### Compilation/Installation
```bash
cd "$HOME/foam/pascal-3.2"
#
./Allwmake.firstInstall
```
