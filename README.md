# Parallel-Programming-Sudoku
Run
To run the solver, supply the executable with command-line arguments for input sudoku file, and modes. I have all the input files under ‘input’ folder.
It supports four modes:
1.	Sequential mode with backtracking algorithm.
2.	Sequential mode with brute force algorithm.
3.	Parallel mode utilizing multithreading.
4.	All modes combined, providing solutions for comparison.
Example of command line arguments: 
./sudoku_file <PATH_TO_INPUT_FILE> <MODE>
./sudoku_file input/16x16_hard.txt 4
Example of how I passed arguments in Visual Studio: ‘Select in VS - Project -> Properties’

 ![alt text](https://github.com/AnaghaDhekne/Parallel-Programming-Sudoku/blob/02edf4b100f39ce7dec0d025c77087cee8e89210/img/Picture1.png)

Output:
 
 
 
 
