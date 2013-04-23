Problem 1
Written by Sondre Engebråten

for (int i = 0; i < n; i++){
	valueToInsert = array[i];
	holePos = i;
	
	while(holePos > 0 && valueToInsert < array[holePos-1]){
		array[holePos] = array[holePos-1];
		holePos = holePos - 1;
	}
	array[holePos] = valueToInsert;

}

Problem 1a

	t0 = 0
startFor:	
	if t0 >= n goto exitFor
	a0 = i*4
	valueToInsert = array[a0];
	holePos = i;
startWhile:	
	t2 = holePos > 0
	a1 = holePos - 1
	a2 = a1*4
	t3 = valueToInsert < array[a2]
	t4 = t2 && t3
	if not t4: goto exitWhile
	a5 = holePos*4
	a6 = holePos - 1
	a7 = a6*4
	array[a5] = array[a7];
	holePos = holePos - 1;
	goto startWhile
exitWhile:
	a8 = holePos * 4
	array[a8] = valueToInsert;
	t0 = t0 + 1
	goto startFor
exitFor:
//}

Problem 1c 


	t0 = 0
startFor:	
	if t0 >= n goto exitFor
	a0 = i*4
	valueToInsert = array[a0];
	holePos = i;
startWhile:	
	t2 = holePos > 0
	a5 = holePos*4
	a2 = a5 - 4
	t3 = valueToInsert < array[a2]
	t4 = t2 && t3
	if not t4: goto exitWhile
	array[a5] = array[a2];
	holePos = holePos - 1;
	goto startWhile
exitWhile:
	array[a5] = valueToInsert;
	t0 = t0 + 1
	goto startFor
exitFor:
//}


Case 1 (Common subexpression):
a5 = holePos*4
a6 = holePos -1
a7 = a6*4
=>
a7 = 4*holePos - 4
=>
a7 = a5 - 4

Case 2 (Dead code):
a6 = holePos - 1
(a6 is never used and can be removed)

Case 3 (Copy propagation):
array[a5] = array[a7] 
=>
array[a5] = array[a2]

Case 4 (Dead code):
a7 = a5 - 4
(a7 is never used)

Case 5 (Common subexpression):
a1 = holePos - 1
a5 = holePos*4
a2 = a1*4
=>
a5 = holePos*4
a2 = a5 - 1
(a1 is only used for a2, by using a5 instead it saves one line)

Case 6 (Copy propagation):
holePos is updated at start of for loop and end of whileloop, however when exiting the loop updating holepos is skipped. As such exitWhile block can reuse holePos from inside the whileloop

a5 = holePos*4
a8 = holePos*4
=>
replace a8 with a5
