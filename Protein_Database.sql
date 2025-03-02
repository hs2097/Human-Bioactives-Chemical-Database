/*Create a database to store Protein structure information*/
CREATE TABLE IF NOT EXISTS PROTEIN_DATABASE(
molname VARCHAR(200),
structure BLOB,
smiles VARCHAR(200),
formula VARCHAR(200),
molwt DECIMAL,
logP DECIMAL,
logD DECIMAL,
numberofatoms INT,
rotatablebonds DECIMAL,
ringcount INT,
H_donors INT,
H_acceptors INT,
valenceelectrons INT,
polarSA DECIMAL,
FAR INT
);





