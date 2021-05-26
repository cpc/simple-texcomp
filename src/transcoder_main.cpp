#include <filesystem>
#include <string>
#include <vector>

#include "transcoder.hpp"

namespace fs = std::filesystem;

/* Exit with error, optionally printing usage */
void err_exit(const std::string& err_msg, bool print_usage=false)
{
    fprintf(stderr, "ERROR: %s\n", err_msg.data());
    if (print_usage)
    {
        fprintf(stderr,
            "Usage: transcoder input_image [<input_image> ...] output_folder\n"
        );
    }
    exit(1);
}

int main(int argc, char **argv)
{
    // Check command line arguments
    if (argc < 3)
    {
        err_exit("Provide at least one image and an output folder", true);
    }

    // Last argument is the output directory
    std::string out_dir = argv[argc-1];
    if ( !(fs::is_directory(out_dir) && fs::exists(out_dir)) )
    {
        err_exit("Can't open output directory", true);
    }

    std::vector<std::string> inp_images(argc-2);
    for (int i = 0; i < argc-2; ++i)
    {
        inp_images.at(i) = std::string(argv[i+1]);
    }

    int ret = simple::transcoder_entry(inp_images, out_dir);
    if (ret != 0)
    {
        err_exit("Failed transcoding files");
    }
}
