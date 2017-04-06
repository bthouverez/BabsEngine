#include "shader.h"

void Shader::destroy()
{
    OpenGLFunctionCore *gl = QOpenGLContext::currentContext()->versionFunctions<OpenGLFunctionCore>();
	gl->glDeleteShader(vertex_);
	gl->glDeleteShader(fragment_);
	gl->glDeleteProgram(program_);
	vertex_ = 0;
	fragment_ = 0;
	program_ = 0;
	initOK_ = false;
}

bool Shader::init()
{
    OpenGLFunctionCore *gl = QOpenGLContext::currentContext()->versionFunctions<OpenGLFunctionCore>();

	//------- Instanciation -------
	if (initShaderStage(vertex_, GL_VERTEX_SHADER, vertexSourceFilename_) == false)
		return false;
	if (initShaderStage(fragment_, GL_FRAGMENT_SHADER, fragmentSourceFilename_) == false)
		return false;

	program_ = gl->glCreateProgram();
	gl->glAttachShader(program_, vertex_);
	gl->glAttachShader(program_, fragment_);
	gl->glLinkProgram(program_);

	//------- Link & check -------
	GLint link(0);
	gl->glGetProgramiv(program_, GL_LINK_STATUS, &link);
	if (link != GL_TRUE)
	{
		GLint err_size(0);
		char *err_msg(NULL);
		gl->glGetProgramiv(program_, GL_INFO_LOG_LENGTH, &err_size);
		err_msg = new char[err_size + 1];
		gl->glGetProgramInfoLog(program_, err_size, &err_size, err_msg);
		err_msg[err_size] = '\0';
		err_out << "--------------------------------------------------------------" << std::endl;
		err_out << "+ GLSL ERROR - linking program + \n" << err_msg << std::endl;
		err_out << "--------------------------------------------------------------" << std::endl;
		delete[] err_msg;
		return false;
	}
	return (initOK_ = true);
}

bool Shader::initShaderStage(GLuint &shader, GLenum shader_stage_type, std::string const &source_filename)
{
    OpenGLFunctionCore *gl = QOpenGLContext::currentContext()->versionFunctions<OpenGLFunctionCore>();

	if (shader_stage_type != GL_VERTEX_SHADER && shader_stage_type != GL_FRAGMENT_SHADER) {
		err_out << "GFX ERROR, Shader stage not supported : " << shader_stage_type << std::endl;
		gl->glDeleteShader(shader);
		return false;
	}
	shader = gl->glCreateShader(shader_stage_type);

	//------- Parse file -------
	std::string code, code_line;
	std::ifstream source_file(source_filename.c_str());
	if (!source_file) {
		err_out << "Error - unable to open file \'" << source_filename << "\' - in gfx::Shader::initShaderStage" << std::endl;
		gl->glDeleteShader(shader);
		return false;
	}
	while (getline(source_file, code_line)) {
		code += code_line + '\n';
	}
	source_file.close();

	//------- Build & check -------
	GLint build_err(0), err_size(0);
	const GLchar* code_cstr = code.c_str();
	gl->glShaderSource(shader, 1, &code_cstr, NULL);
	gl->glCompileShader(shader);
	gl->glGetShaderiv(shader, GL_COMPILE_STATUS, &build_err);
	if (build_err != GL_TRUE)
	{
		gl->glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &err_size);
		char *err_msg = new char[err_size + 1];
		err_msg[err_size] = '\0';
		gl->glGetShaderInfoLog(shader, err_size, &err_size, err_msg);
		err_out << "--------------------------------------------------------------" << std::endl;
		err_out << "+ GLSL ERROR - shader stage (" << source_filename << "::" << shader_stage_type << ") +" << std::endl;
		err_out << "--------------------------------------------------------------" << std::endl;
		err_out << err_msg << std::endl;
		err_out << "--------------------------------------------------------------" << std::endl;
		delete[] err_msg;
		return false;
	}

	//--- OK ---
	return true;
}

Shader& Shader::operator=(Shader const &shader)
{
	if (this != &shader) {
		vertexSourceFilename_ = shader.vertexSourceFilename_;
		fragmentSourceFilename_ = shader.fragmentSourceFilename_;
		initOK_ = shader.initOK_;
		destroy();
		init();
	}
	return *this;
}

GLuint Shader::getProgramID() const {
	return program_;
}