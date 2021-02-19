# Icosahedral matices generator for protein crystallographer

## preparation

After getting the file, make it executable.
```
chmod 755 ./icosagenmat.py
```

## how to run

```
./icosagenmat.py 

./icosagenmat.py 90 0 0 

./icosagenmat.py 20 30 10 115.2 123.1 134.1

./icosagenmat.py 20 30 10 115.2 123.1 134.1 "FIXR ON FIXB ON" 
```

## explanation
- Matrices numbering is the same as viperdb.scripps.edu one 
- inputs are alpha beta gamma x_center y_center z_center in Euler angles and in orthogonal coordinate.
- Units for input angle angles are degree, coordinates for center are in orthogonal as x y z 
- Outputs are phaser_sol_eul refmac_ncscon_eul  refmac_ncscon_mat and biomat
- After 6 numbers of angles and coordinates, the characters with " " can be inputted for phaser input file.


