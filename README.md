# Voronoi Maze

Generate Voronoi Maze to out.txt file by running:
```
g++ .\mazeGen.cpp .\efficientPoisson.cpp
./a.exe > out.txt
```
Format of out.txt:
1) First two lines are coordinates of the start room and the end room.
2) The remaining lines are the coordinate of the walls

To visualize the maze:
1) Create Virtual Env by:
   ```
   python -m venv env
   ```
2) Install requirements:
   ```
   .\env\Scripts\activate (For windows)
   pip install -r requirements.txt
   ```
3) Simulate the Maze:
   ```
   python .\voronoi.py
   ```
