# cpp_phong_rendering
Using cpp to perform what is called Blinn-Phong Shading rendering to the scene which is created totally using cpp ...... The scene specs including (shapes, lightsources, backgroundcolor, depth...etc) is read from a JSON file and then the rendering function works to render the scene based on the algorithm of blinn-phong .... 

<p style="text-align: center;"> Blinn-Phong rendering final image with texture </p> <br>
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/24fda539-20fd-447c-99b6-ffbea22e935f" alt="phong_scene_textured" width="400" style="text-align: center;">


## Project process
To know more about phong rendering visit this website... [Phong Shading](https://en.wikipedia.org/wiki/Phong_shading)
Understanding the concept of phong rendering well and the structure of the code, I managed to use some of the LLM sources, like ChatGPT and Copilot, to fully create the whole rendering project while adding my own specs to the rendering process, like the colored binary rendering. 

## specs added to the rendering process (Feature list):
- Image write (ppm file) without using any helping libraries in cpp project I extract the images in a ppm format, Portable Pixmap.
- Pinhole camera implementation
- Intersection between different shapes (sphere, cylinder, triangle)
- Binary rendering
- Shadows
- Textures are added to the material (wooden texture)
- Tone mapping
- Reflections
- Refraction
- Multithreading
- Multisampling per pixel
- Defocus in finite aperture cameras
- Multi-bounce path tracing
- A video showcasing the quality of rendering process.

## compiler 
if you are using gcc compiler then make sure to go to the directory of the src code which is ./code and run the following command after adding the json library bin folder to your environment variables. 
<br> ` g++ main.cc -o main -I .\json-develop\json-develop\include. `
The current program runs blinn-phong rendering process for one of the json files called "scene.json."
More options to run binary rendering on the same file, "scene.json," or even to different scenes' files.

<br> and then you will be able to run the .exe file main.exe

 ### Linux users 
 For Linus users, there is a Makefile to use so, simply you can open the command line in the ./code directory and run the command
 `make`

 ## video for the final phong rendering process 
https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/2b20cfa9-666a-4955-a921-e756b4f52064



 ## images for phong rendering process 

 final phong rendering image with added wooden texture to the ground ....
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/24fda539-20fd-447c-99b6-ffbea22e935f" alt="phong_scene_textured" width="400">


Intermediary step for phong rendering before edits and without texture... 
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/d62dd3ec-f6e6-4962-9df0-914a8d2cd298" alt="phong_scene_textured" width="400" style="text-align: center;">


phong rendering for another scene. 
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/60ece903-b33f-41a8-b891-369e8e79a1ee" alt="phong_scene_textured" width="400" style="text-align: center;">


## images for binary rendering process

final binary rendering images for different scenes and json files in the project <br>
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/9c244922-b4dc-41d8-b1aa-64a3256bca83" alt="phong_scene_textured" width="400" style="text-align: center;">
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/83bbd5ba-88a0-482d-b6e5-21116b100daa" alt="phong_scene_textured" width="400" style="text-align: center;">
<img src="https://github.com/moneerzaki/Rendering_Modeling_cpp/assets/78418503/78008341-3eff-479f-9282-93dd2cbcaf79" alt="phong_scene_textured" width="400" style="text-align: center;">




