# For Linux :

## How to install (the whole project) :

```
git clone -b invariant https://github.com/soufianeelm/PetriSpot.git ;
cd PetriSpot ;
./build.sh
```


## How to install (just the program) :

LA PREMIERE COMMANDE NE RAMENE PAS LE BON FILE

```
wget --progress=dot:mega https://github.com/soufianeelm/PetriSpot/blob/Inv-Linux/petri ;
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

--Psemiflow || --Pflow

--Tsemiflow || --Tflow
