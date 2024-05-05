# For Linux :

## How to install (the whole project) :

```
git clone -b invariant https://github.com/soufianeelm/PetriSpot.git ;
cd PetriSpot ;
./build.sh ;
cd Petri/src
```


## How to install (just the program) :

```
git clone -b Inv-Linux https://github.com/soufianeelm/PetriSpot.git ;
cd PetriSpot ;
```


## How to use :

```
./petri -i [model.pnml] [flags]
```


### [model.pnml] : 

path of model file in format .pnml


### [flags] :

-q (quiet)

--findDeadlock

--Psemiflows || --Pflows

--Tsemiflows || --Tflows
