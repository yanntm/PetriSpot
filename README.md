# For Linux :

## How to install (the whole project) :

```
git clone -b invariant https://github.com/soufianeelm/PetriSpot.git

cd PetriSpot

./build.sh
```


## How to install (just the program) :

```
wget --progress=dot:mega https://github.com/soufianeelm/PetriSpot/blob/Inv-Linux/petri

chmod +x petri
```


## How to use :

```
cd Petri/src

./petri [model.pnml] [flags]
```


### [model.pnml] : 

path of model file in format .pnml


### [flags] :

--findDeadlock

--Psemiflows || --Pflows

--Tsemiflows || --Tflows
