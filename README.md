# cpp_phong_rendering
using cpp to perform what is called blinn-phong rendering to the scene which is created totally using cpp ...... The scene specs including (shapes, lightsources, backgroundcolor, depth...etc) is read from a JSON file and then the rendering function works to render the scene based on the algorithm of blinn-phong .... 

# compiler 
if you are using gcc compiler then make sure to go to the directory of the src code which is ./code and run the following command after adding the json library bin folder to your environment variables. 
<br> ` g++ main.cc -o main -I .\json-develop\json-develop\include. `
The current program runs blinn-phong rendering process for one of the json files called "scene.json."
More options to run binary rendering on the same file, "scene.json," or even to different scenes' files.

<br> and then you will be able to run the .exe file main.exe

 # Linux users 
 For Linus users, there is a Makefile to use so, simply you can open the command line in the ./code directory and run the command
 `make`

 # images for phong rendering process 
 

