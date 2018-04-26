BOOST_INCLUDE=~/boost_1_57_0/boost

.PHONY: all random_walk random_walk_wikipedia

all: random_walk random_walk_wikipedia

random_walk:
	g++ -g -std=c++11 -O3 -fopenmp -I$(BOOST_INCLUDE) -I include -lboost_random -ltcmalloc_minimal random_walk.cpp -o random_walk

random_walk_wikipedia:
	g++ -g -std=c++11 -O3 -fopenmp -I$(BOOST_INCLUDE) -I include -lboost_random -ltcmalloc_minimal random_walk_wikipedia.cpp -o random_walk_wikipedia

.PHONY: clean
clean:	
	rm random_walk random_walk_wikipedia
