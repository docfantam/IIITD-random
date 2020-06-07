ks: ks.cpp data.txt
	g++ ks.cpp -o ks -lgmp

data.txt: create_data
	./create_data

create_data: create_data.cpp
	g++ create_data.cpp -o create_data -lgmp
