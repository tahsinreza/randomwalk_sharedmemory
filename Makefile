BOOST_INCLUDE=~/boost_1_57_0/boost

.PHONY: random_walk clean

random_walk:
	g++ -g -std=c++11 -O3 -fopenmp -I$(BOOST_INCLUDE) -I include -lboost_random -ltcmalloc_minimal random_walk.cpp -o random_walk

clean:	
	rm random_walk
