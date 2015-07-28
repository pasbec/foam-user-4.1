# README
------

## Installation
------

### Pre-requirements

#### Debian
##### Debian Wheezy ([Ubuntu 12.04](#ubuntu1204))
##### Debian Sid ([Ubuntu 13.10](#ubuntu1310))
##### Debian Jessie ([Ubuntu 14.04](#ubuntu1404))

#### Fedora
##### Fedora 19
```bash
sudo yum groupinstall "Development Tools"
sudo yum install gcc-c++ binutils-devel bison flex m4 zlib-devel \
qt-devel qtwebkit-devel mercurial graphviz
```
##### Fedora 20
```bash
sudo yum groupinstall "Development Tools"
sudo yum install gcc-c++ binutils-devel cmake wget bison flex m4 zlib-devel \
qt-devel qtwebkit-devel mercurial graphviz
```

#### openSUSE
##### openSUSE 13.1
```bash
su -
zypper install git openmpi-devel gnuplot -t pattern devel_qt4
exit
```

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
git subtree pull -P 'pascal-3.1/etc/bashrc.d' 'dotfiles' 'subtree_dot-bashrc.d' --squash
```
#### Update submodule foam-extend
```bash
cd "$HOME/foam"
#
git submodule update --remote 'foam-extend-3.1'
```
#### User data for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
git config user.name 'Pascal Beckstein'
git config user.email 'pascal@becksteins.de'
```
#### Git repositories for submodule foam-extend (optional)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
git remote add 'upstream' 'https://git.code.sf.net/p/openfoam-extend/foam-extend-3.1'
git remote add 'upstream_mirror' 'https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-foam-extend-3.1.git'
git remote add 'sf' 'ssh://p-be@git.code.sf.net/u/p-be/openfoam-extend'
#git remote set-url 'origin' 'https://git.becksteins.de/foam-extend/foam-extend-3.1'
```


### Project part

#### Official documentation
```bash
cd "$HOME/foam/foam-extend-3.1"
#
less 'doc/buildInstructions/Ubuntu/Ubuntu_14.04'
```
#### Source environment
```bash
cd "$HOME/foam/foam-extend-3.1"
#
source 'etc/bashrc'
```
#### Prepare ParaView compilation (qmake)
```bash
cd "$HOME/foam/foam-extend-3.1"
#
export QT_BIN_DIR="$(dirname $(which qmake))"
echo "export QT_BIN_DIR=$(dirname "$(which qmake)")" >> 'etc/prefs.sh'
cd -
```
#### Compilation/Installation
```bash
cd "$HOME/foam/foam-extend-3.1"
#
./Allwmake.firstInstall
```


### User part

#### Source environment
```bash
cd "$HOME/foam/pascal-3.1"
#
source 'etc/bashrc'
fe31
```
#### Compilation/Installation
```bash
cd "$HOME/foam/pascal-3.1"
#
./Allwmake
```