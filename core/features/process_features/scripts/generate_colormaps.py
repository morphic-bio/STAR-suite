import numpy as np
import matplotlib.pyplot as plt

def generate_colormap(name, num_colors):
    """
    Generate a colormap with the specified name and number of colors.

    Args:
        name (str): Name of the colormap (e.g., "viridis", "plasma").
        num_colors (int): Number of colors to sample.

    Returns:
        np.ndarray: Array of RGB colors.
    """
    cmap = plt.cm.get_cmap(name, num_colors)  # Get the colormap
    colors = cmap(np.linspace(0, 1, num_colors))[:, :3]  # Extract RGB values
    return colors

def format_as_c_array(colors, colormap_name):
    """
    Format the colormap as a C array.

    Args:
        colors (np.ndarray): Array of RGB colors.
        colormap_name (str): Name of the colormap.

    Returns:
        str: C code defining the colormap as a constant array.
    """
    array_name = f"{colormap_name}_colormap_{len(colors)}"
    c_code = f"const double {array_name}[{len(colors)}][3] = {{\n"
    for color in colors:
        r, g, b = color
        c_code += f"    {{{r:.6f}, {g:.6f}, {b:.6f}}},\n"
    c_code = c_code.rstrip(",\n") + "\n};\n"  # Remove the trailing comma
    return c_code

def save_colormap_to_file(colors, colormap_name, filename):
    """
    Save the colormap as a C array to a file.

    Args:
        colors (np.ndarray): Array of RGB colors.
        colormap_name (str): Name of the colormap.
        filename (str): Output file name.
    """
    c_code = format_as_c_array(colors, colormap_name)
    with open(filename, "w") as file:
        file.write(c_code)

if __name__ == "__main__":
    # Define resolutions and colormap name
    resolutions = [16, 64, 256, 1024]
    colormap_name = "plasma"  # Change to "viridis" or other colormaps if needed

    for res in resolutions:
        # Generate the colormap
        colors = generate_colormap(colormap_name, res)

        # Save to a file
        filename = f"{colormap_name}_colormap_{res}.h"
        save_colormap_to_file(colors, colormap_name, filename)
        print(f"Saved: {filename}")
