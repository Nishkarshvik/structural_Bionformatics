                                                    
                                                    #################################################################################################
                                                    #                                                                                               #                                                                                                           #                       Python script for converting SVG To HIGH RESOLUTION PNG                 #                                                                                                           #                                                                                               #                                                                                                           #                          DEVELOPED BY: Niskarsh Vikram Singh                                  #                                                                                                           #                                JRF at CARI, Kolkata                                           #                                                                                                           #                                                                                               #                                                                                                           #                            Cytoscape SVG Image conversion                                     #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #                                                                                               #                                                                                                           #################################################################################################                                                                                                           


# script can be used through the command line by providing file_name.svg with dpi as per the desire.
# usage: script_name file_name.svg dpi
# https://github.com/Nishkarshvik/structural_Bionformatics


import sys
import os
import cairosvg

def convert_svg_to_png(input_file, dpi):
    """
    Converts an SVG file to a high-resolution PNG.

    Args:
        input_file (str): Path to the SVG file.
        dpi (int): Desired resolution in dots per inch (DPI).

    Returns:
        str: Path to the generated PNG file.
    """
    try:
        # Check if the file exists
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"The file '{input_file}' does not exist.")

        # Ensure the file has an .svg extension
        if not input_file.lower().endswith('.svg'):
            raise ValueError("The input file must have an .svg extension.")

        # Generate output file name
        output_file = os.path.splitext(input_file)[0] + ".png"

        # Read the SVG file dimensions
        with open(input_file, "r") as file:
            svg_content = file.read()
        width = int(svg_content.split("width=\"")[1].split("\"")[0].replace("px", ""))
        height = int(svg_content.split("height=\"")[1].split("\"")[0].replace("px", ""))

        # Calculate scaling factor based on DPI
        scale_factor = dpi / 96  # 96 DPI is the standard resolution for SVGs
        new_width = int(width * scale_factor)
        new_height = int(height * scale_factor)

        # Convert SVG to PNG with specified width and height
        cairosvg.svg2png(url=input_file, write_to=output_file, output_width=new_width, output_height=new_height)
        print(f"Conversion successful! High-resolution PNG saved as: {output_file}")
        return output_file

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    import argparse

    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Convert an SVG file to a high-resolution PNG image.")
    parser.add_argument("input_file", type=str, help="Path to the input SVG file.")
    parser.add_argument("dpi", type=int, help="Desired resolution in DPI (dots per inch).")

    # Parse the arguments
    args = parser.parse_args()

    # Convert the SVG to PNG
    convert_svg_to_png(args.input_file, args.dpi)
