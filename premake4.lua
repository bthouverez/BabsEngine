-- A solution contains projects, and defines the available configurations
solution "RenderBoy"

   root_dir = path.getabsolute("./")
   src_dir  = path.getabsolute("./src/")
   tests_dir = path.getabsolute("./tests/")

   includedirs { ".", "./src/" }

   configurations { "Debug", "Release" }
 
   -- A project defines one build target
   project "RenderBoy"
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
tests = {
   "helloworld_c++14",
   "helloworld_c++11"
}


for i, name in ipairs(tests) do
   project(name)
      language "C++"
      kind "ConsoleApp"
      if name == "helloworld_c++14" then
         buildoptions { "-std=c++14" }
      elseif name == "helloworld_c++11" then
         buildoptions { "-std=c++11" }
      end
      targetdir "bin"
      files { tests_dir .. "/" .. name .. "/main.cpp" }
end