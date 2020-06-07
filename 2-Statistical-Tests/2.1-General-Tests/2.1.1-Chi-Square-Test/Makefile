cs: cs.cpp data.txt
	g++ cs.cpp -o cs -lgmp

data.txt: create_data
	./create_data

create_data: create_data.cpp
	g++ create_data.cpp -o create_data -lgmp
