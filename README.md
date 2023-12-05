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

 final phong rendering image with added wooden texture to the ground ....
![phong_scene_textured](https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/24fda539-20fd-447c-99b6-ffbea22e935f)

Intermediary step for phong rendering before edits and without texture... 
![phong_scene_before](https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/d62dd3ec-f6e6-4962-9df0-914a8d2cd298)

phong rendering for another scene. 
![simple_phong_rendering](https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/60ece903-b33f-41a8-b891-369e8e79a1ee)

# images for binary rendering process

final binary rendering images for different scenes and json files in the project 
![binary_mirror_image](https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/9c244922-b4dc-41d8-b1aa-64a3256bca83)
![binary_scene](https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/83bbd5ba-88a0-482d-b6e5-21116b100daa)
![binary_simple_phong](https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/78008341-3eff-479f-9282-93dd2cbcaf79)



