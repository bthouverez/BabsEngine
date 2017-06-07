-- A solution contains projects, and defines the available configurations
solution "RenderBoy"

   root_dir = path.getabsolute("./")
   src_dir  = path.getabsolute("./src/")
   tests_dir = path.getabsolute("./tests/")

   includedirs { ".", "./src/" }

   configurations { "Debug", "Release" }
 
   -- A project defines one build target
   project "RenderBoy"
      location "build"
      kind "ConsoleApp"
      language "C++"
      files { "src/**.h", "src/**.cpp" }
 
      configuration "Debug"
         defines { "DEBUG" }
         flags { "Symbols" }
         buildoptions { "-std=c++14" }
         links { "GLEW", "SDL2", "GL" }
         targetdir "bin"
 
      configuration "Release"
         defines { "NDEBUG" }
         flags { "Optimize" }
         buildoptions { "-std=c++14" }
         links { "GLEW", "SDL2", "GL" }
         targetdir "bin"

-- Tests
-- To create a test, create a directory in 'tests' directory and create a main.cpp file into
-- Then add the name of the directory you just created in the following list
tests = {
   "helloworld_c++14",
   "helloworld_c++11",
   "test_vec3"
}

--[[
project("helloworld_c++14")
   language "C++"
   kind "ConsoleApp"
   targetdir "bin"
   location "build-tests"
   buildoptions { "-std=c++14" }
   files { tests_dir .. "/helloworld_c++14/main.cpp" }

project("helloworld_c++11")
   language "C++"
   kind "ConsoleApp"
   targetdir "bin"
   location "build-tests"
   buildoptions { "-std=c++11" }
   files { tests_dir .. "/helloworld_c++11/main.cpp" }

project("test_vec3")
   language "C++"
   kind "ConsoleApp"
   targetdir "bin"
   location "build-tests"
   files { tests_dir .. "/test_vec3/main.cpp", src_dir .. "/Algebra/vec3.h", src_dir .. "/Algebra/vec3.cpp"}
--]]

for i, name in ipairs(tests) do
   project(name)
      language "C++"
      kind "ConsoleApp"
      if name == "helloworld_c++14" then
         buildoptions { "-std=c++14" }
      elseif name == "helloworld_c++11" then
         buildoptions { "-std=c++11" }
      elseif name == "test_vec3" then
         buildoptions { "-std=c++11" }
         files { src_dir .. "/Algebra/vec3.h", src_dir .. "/Algebra/vec3.cpp", }
      end
      targetdir "bin"
      location "build-tests"
      files { tests_dir .. "/" .. name .. "/main.cpp" }
end
