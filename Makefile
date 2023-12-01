FILE = expected_iter

all:
	g++ -Wall -Wextra -g -fsanitize=address,undefined -D_GLIBCXX_DEBUG -o $(FILE) $(FILE).cpp

fast:
	g++ -O3 -march=native -o $(FILE) $(FILE).cpp