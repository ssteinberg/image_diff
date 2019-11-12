
/*
*
* Image diff tool
* Copyright  Shlomi Steinberg
*
*/

#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <cassert>
#include <getopt.h>
#include <vector>

#include <png.h>

// color-util expects M_PI to be defined, for reasons
#define M_PI glm::pi<double>()
#include <color-util/XYZ_to_Lab.hpp>
#include <color-util/CIEDE2000.hpp>
#undef M_PI

#include <glm/glm.hpp>

struct surface {
	std::unique_ptr<std::uint8_t[]> data;
	std::size_t width{}, height{}, components{};
};

surface load_png(const std::filesystem::path &file_name) {
	png_byte header[8];

	FILE *fp = fopen(file_name.string().data(), "rb");
	if (!fp) {
		throw std::runtime_error("Could not open file");
	}

	// read the header
	fread(header, 1, 8, fp);

	if (png_sig_cmp(header, 0, 8)) {
		std::cerr << file_name << " is not a PNG" << std::endl;
		fclose(fp);
		throw std::runtime_error("Not a valid PNG");
	}

	png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png_ptr) {
		std::cerr << file_name << " png_create_read_struct returned 0" << std::endl;
		fclose(fp);
		throw std::runtime_error("Not a valid PNG");
	}

	// create png info struct
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		std::cerr << file_name << " png_create_info_struct returned 0" << std::endl;
		png_destroy_read_struct(&png_ptr, nullptr, nullptr);
		fclose(fp);
		throw std::runtime_error("Not a valid PNG");
	}

	// create png info struct
	png_infop end_info = png_create_info_struct(png_ptr);
	if (!end_info) {
		std::cerr << file_name << " png_create_info_struct returned 0" << std::endl;
		png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
		fclose(fp);
		throw std::runtime_error("Not a valid PNG");
	}

	// the code in this if statement gets called if libpng encounters an error
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr << file_name << " error from libpng" << std::endl;
		png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
		fclose(fp);
		throw std::runtime_error("libpng error");
	}

	// init png reading
	png_init_io(png_ptr, fp);

	// let libpng know you already read the first 8 bytes
	png_set_sig_bytes(png_ptr, 8);

	// read all the info up to the image data
	png_read_info(png_ptr, info_ptr);

	// variables to pass to get info
	int bit_depth, color_type;
	png_uint_32 temp_width, temp_height;

	// get info about png
	png_get_IHDR(png_ptr, info_ptr, &temp_width, &temp_height, &bit_depth, &color_type,
				 nullptr, nullptr, nullptr);

	if (bit_depth != 8 && (bit_depth != 1 || color_type != PNG_COLOR_TYPE_GRAY)) {
		std::cerr << file_name << " Unsupported bit depth " << std::to_string(bit_depth) << ".  Must be 8" << std::endl;
		throw std::runtime_error("Unsupported bit depth");
	}

	int components;
	switch (color_type) {
	case PNG_COLOR_TYPE_GRAY:
		components = 1;
		break;
	case PNG_COLOR_TYPE_RGB:
		components = 3;
		break;
	case PNG_COLOR_TYPE_RGB_ALPHA:
		components = 4;
		break;
	default:
		std::cerr << file_name << " Unknown libpng color type " << std::to_string(color_type) << std::endl;
		fclose(fp);
		throw std::runtime_error("Unsupported PNG color type");
	}

	// Update the png info struct.
	png_read_update_info(png_ptr, info_ptr);

	// Row size in bytes.
	auto rowbytes = png_get_rowbytes(png_ptr, info_ptr);
	// Allocate the image_data as a big block
	const size_t w = rowbytes / components;
	const auto level0_size = rowbytes * temp_height;
	std::unique_ptr<uint8_t[]> image_data = std::make_unique<uint8_t[]>(level0_size);

	// row_pointers is for pointing to image_data for reading the png with libpng
	png_byte ** row_pointers = reinterpret_cast<png_byte **>(malloc(temp_height * sizeof(png_byte *)));
	if (bit_depth == 8) {
		// set the individual row_pointers to point at the correct offsets of image_data
		for (unsigned int i = 0; i < temp_height; i++)
			row_pointers[i] = image_data.get() + i * w * components;

		// read the png into image_data through row_pointers
		png_read_image(png_ptr, row_pointers);
	}
	else {
		std::cerr << file_name << " Unsupported bit depth" << std::endl;
		fclose(fp);
		throw std::runtime_error("Unsupported PNG bit depth");
	}

	fclose(fp);
	png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
	free(row_pointers);

	return surface{ std::move(image_data), w, temp_height, static_cast<std::size_t>(components) };
}

void write_2d(const char *file_name, const surface &s) {
	const uint8_t *image_data = s.data.get();
	int components = s.components, width = static_cast<int>(s.width), height = static_cast<int>(s.height);

	FILE *fp = fopen(file_name, "wb");
	if (!fp) {
		std::cerr << file_name << " can't be opened for writing" << std::endl;
		throw std::runtime_error("PNG save error");
	}

	png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
	if (!png) {
		std::cerr << file_name << " png_create_write_struct failed" << std::endl;
		fclose(fp);
		throw std::runtime_error("PNG save error");
	}

	png_infop info = png_create_info_struct(png);
	if (!info) {
		std::cerr << file_name << " png_create_info_struct failed" << std::endl;
		fclose(fp);
		throw std::runtime_error("PNG save error");
	}

	if (setjmp(png_jmpbuf(png))) {
		std::cerr << file_name << " png_jmpbuf failed" << std::endl;
		fclose(fp);
		throw std::runtime_error("PNG save error");
	}

	png_byte ** const row_pointers = reinterpret_cast<png_byte **>(malloc(height * sizeof(png_byte *)));
	if (row_pointers == nullptr) {
		std::cerr << file_name << " could not allocate memory for PNG row pointers" << std::endl;
		fclose(fp);
		throw std::runtime_error("PNG save error");
	}

	// set the individual row_pointers to point at the correct offsets of image_data
	// To maintain compatibility png_write_image requests a non-const double pointer, hack the const away...
	for (int i = 0; i < height; i++)
		row_pointers[i] = const_cast<png_byte*>(image_data + i * width * components);

	png_init_io(png, fp);

	int color_type;
	switch (components) {
	case 1: color_type = PNG_COLOR_TYPE_GRAY; break;
	case 3: color_type = PNG_COLOR_TYPE_RGB; break;
	case 4: color_type = PNG_COLOR_TYPE_RGBA; break;
	default: {
		std::cerr << "Invalid component count" << std::endl;
		throw std::runtime_error("PNG save error");
	}
	}
	png_set_IHDR(png,
		info,
		width, height,
		8,
		color_type,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);
	png_write_info(png, info);

	png_write_image(png, row_pointers);
	png_write_end(png, nullptr);

	free(row_pointers);
	fclose(fp);
}

std::vector<double> readCSVline(std::istream& str) {
	std::vector<double> result;
	std::string l;
	std::getline(str,l);

	std::stringstream ss(l);
	std::string s;
	while (std::getline(ss, s, ','))
		result.push_back(std::stod(s));
	return result;
}

float de00_error_color_code(float de00) {
	const float A = 0.75f;
	const float B = 1.f;
	const float M = 8.f;

	if (de00 <= A)
		return de00 * .05f;
	else if (de00 <= B) {
		const auto m = (de00 - A) / (B - A);
		return glm::log(glm::mix(1.f, 2.f, m))*0.3f + de00_error_color_code(A);
	}
	else
		return glm::pow((de00 - B) / (M - B), 2.3f) + de00_error_color_code(B);
}

void print_usage() {
	std::cout << "Usage: image_diff <options> --out=<output> --map=<map name> <input1> <input2>" << std::endl;
	std::cout << "-m <map>\t\tSpecify colormap." << std::endl;
	std::cout << "-s <scale>\t\tScale errors." << std::endl;
	std::cout << "-R\t\t\t\tAssume inputs are in RGB space, as opposed to sRGB." << std::endl;
}

int main(int argc, char* argv[]) {
	const char* path = nullptr, *map = "cividis.csv", *input1 = nullptr, *input2 = nullptr;
	bool srgb = true;
	float scale = 1.f;

	if (argc > 3) {
		input1 = argv[argc-2];
		input2 = argv[argc-1];

		option long_options[] = {
			{"out", required_argument, NULL, 'o'},
			{"map", required_argument, NULL, 'm'},
			{"scale", required_argument, NULL, 's'},
			{"nosrgb", no_argument, NULL, 'r'},
			{NULL, 0, NULL, 0}
		};
		char ch;
		while ((ch = getopt_long(argc-2, argv, "o:m:s:R", long_options, NULL)) != -1) {
			switch (ch) {
			case 'o':
				path = optarg;
				break;
			case 'm':
				map = optarg;
				break;
			case 's':
				scale = static_cast<float>(atof(optarg));
				break;
			case 'R':
				srgb = false;
				break;
			default:
				std::cerr << "Unkown option!" << std::endl;
				return 1;
			}
		}
	}

	if (!path || !map || !input1 || !input2 || scale<=.0) {
		print_usage();
		return 1;
	}
	
	const std::filesystem::path map_path{ map };
	const std::filesystem::path input_paths[2] = { input1, input2 };

	for (const auto &p : input_paths)
		if (!std::filesystem::exists(p)) {
			std::cerr << "'" << p.string() << "' not found!" << std::endl;
			return 1;
		}

	// Read map
	std::vector<glm::vec3> colormap;
	{
		std::ifstream ifs(map_path);
		if (!ifs) {
			std::cerr << "'" << map_path << "' not found!" << std::endl;
			return 1;
		}
		std::string dummy;
		std::getline(ifs, dummy);
		std::vector<double> data;
		colormap.reserve(300);
		while ((data = readCSVline(ifs)).size())
			colormap.emplace_back(glm::vec3{ static_cast<float>(data[0]), static_cast<float>(data[1]), static_cast<float>(data[2]) });
	}

	// Read PNGs
	surface ins[2];
	for (int i=0;i<2;++i)
		ins[i] = load_png(input_paths[i]);

	// Sanity
	if (!ins[0].width || !ins[0].height || ins[0].width!=ins[1].width || ins[0].height!=ins[1].height) {
		std::cerr << "Inputs must have identical dimensions" << std::endl;
		return 1;
	}

	surface out;
	out.width = ins[0].width;
	out.height = ins[0].height;
	out.components = 3;
	out.data = std::make_unique<std::uint8_t[]>(out.width*out.height*out.components);

	// Calculate diff
	const auto RGBtoXYZ = glm::inverse(glm::mat3( 2.0413690, -0.5649464, -0.3446944,
												 -0.9692660,  1.8760108,  0.0415560,
												  0.0134474, -0.1183897,  1.0154096));
	const auto pixels = out.width * out.height;
	for (int x=0;x<pixels;++x) {
		auto srcRGB = ins[0].components >= 3 ? 
			glm::vec3{ static_cast<float>(ins[0].data[x*ins[0].components + 0]), static_cast<float>(ins[0].data[x*ins[0].components + 1]), static_cast<float>(ins[0].data[x*ins[0].components + 2]) } / 255.f : 
			glm::vec3{ static_cast<float>(ins[0].data[x]) } / 255.f;
		auto dstRGB = ins[1].components >= 3 ? 
			glm::vec3{ static_cast<float>(ins[1].data[x*ins[1].components + 0]), static_cast<float>(ins[1].data[x*ins[1].components + 1]), static_cast<float>(ins[1].data[x*ins[1].components + 2]) } / 255.f : 
			glm::vec3{ static_cast<float>(ins[1].data[x]) } / 255.f;
		if (srgb) {
			srcRGB = glm::pow(srcRGB, glm::vec3{ 2.2f });
			dstRGB = glm::pow(dstRGB, glm::vec3{ 2.2f });
		}

		const auto srcXYZ = srcRGB * RGBtoXYZ;
		const auto dstXYZ = dstRGB * RGBtoXYZ;
		const auto srclab = colorutil::convert_XYZ_to_Lab(colorutil::XYZ{ srcXYZ.x, srcXYZ.y, srcXYZ.z });
		const auto dstlab = colorutil::convert_XYZ_to_Lab(colorutil::XYZ{ dstXYZ.x, dstXYZ.y, dstXYZ.z });

		const auto ciede2000 = colorutil::calculate_CIEDE2000(srclab, dstlab) * scale;

		const auto cmidxf = de00_error_color_code(ciede2000) * colormap.size();
		const auto cm0 = colormap[std::min(static_cast<size_t>(cmidxf),   colormap.size() - 1)];
		const auto cm1 = colormap[std::min(static_cast<size_t>(cmidxf)+1, colormap.size() - 1)];
		const auto cm = glm::mix(cm0, cm1, glm::fract(cmidxf));

		const auto cm_srgb = glm::pow(cm, glm::vec3{ 1/2.2f });
		const auto packed = glm::packUnorm4x8(glm::vec4{ cm_srgb, .0f });

		out.data[x*3 + 0] = static_cast<uint8_t>((packed >> 0) & 0xFF);
		out.data[x*3 + 1] = static_cast<uint8_t>((packed >> 8) & 0xFF);
		out.data[x*3 + 2] = static_cast<uint8_t>((packed >> 16) & 0xFF);
	}

	// Write result
	write_2d(path, out);

	return 0;
}
