#ifndef SHADER_HEADER
#define SHADER_HEADER

#include <iostream>
#include <string>
#include <fstream>
#include "gl_version.h"

class Yann_Shader
{
public:
    Yann_Shader() : err_out(std::cerr) { }

    Yann_Shader(std::string vertexSourceFilename, std::string fragmentSourceFilename, std::ostream & error_out)
		: err_out(error_out)
		, vertexSourceFilename_(vertexSourceFilename)
		, fragmentSourceFilename_(fragmentSourceFilename)
		, vertex_(0), fragment_(0), program_(0), initOK_(false)
        { std::cout << "Shader construction" << std::endl; }

    Yann_Shader(Yann_Shader const &shader)
		: err_out(std::cerr)
		, vertexSourceFilename_(shader.vertexSourceFilename_)
		, fragmentSourceFilename_(shader.fragmentSourceFilename_)
		, vertex_(0), fragment_(0), program_(0), initOK_(shader.initOK_)
	{
		if (initOK_)
			init();
	}

    ~Yann_Shader()
	{
		destroy();
	}

	bool init();
	void destroy();
	GLuint getProgramID() const;
    Yann_Shader& operator=(Yann_Shader const &shader);


private:
	bool initShaderStage(GLuint &shader, GLenum shader_stage_type, std::string const &source);

private:
	std::ostream & err_out;
	std::string vertexSourceFilename_;
	std::string fragmentSourceFilename_;
	GLuint vertex_;
	GLuint fragment_;
	GLuint program_;
	bool initOK_;
};


#endif // SHADER_HEADER
