import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter


def con_line(x_mesh, y_mesh, z, clim):
        # 5. control the contour line _con_line_.  
    if con_line is None:
        pass 
    else:
        if con_line[2] == 0:
            if clim is None:
                contour_lines = plt.contour(x_mesh, y_mesh, z, colors='black', linewidths=0.5)
            else: 
                contour_levels1 = list(np.linspace(clim[0], 0, 6))
                contour_levels2 = list(np.linspace(0, clim[1], 6))
                contour_levels = contour_levels1 + contour_levels2[1:]
                contour_lines = plt.contour(x_mesh, y_mesh, z, colors='black', levels=contour_levels, linewidths=0.5) # @ contour lines
        else:
            contour_levels = np.linspace(con_line[0], con_line[1], con_line[2])
            contour_lines = plt.contour(x_mesh, y_mesh, z, levels=contour_levels, colors='black', linewidths=0.5) # @ contour lines

# 1. Plot function: contour plot. 
def contour_plot(x, y, z, pltfig=None, label=None, ft=None, MyColor="PiYG", 
                 xylog=None, clim=None, con_line=None, pltcbar=1, ax_line=None, 
                 cticks=None, xticks=None, yticks=None, grid=None, adjust=None,
                 ContourLable=None):
    """Update: 2025.12.30."""

    # 2. Some default settings: figure size, fontsize, and lables. 
    if pltfig is None:
        plt.figure(figsize=(8, 6))
    else:
        if len(pltfig) == 2:
            plt.figure(figsize=(pltfig[0], pltfig[1]))
        else:
            pass


    # 3. Set the label, font size, and font name.
    if label is None:
        xname, yname, title, cbar_name = "x", "y", "Title", "color bar"
    else:
        xname, yname, title, cbar_name = label[0], label[1], label[2], label[3]
    if ft is None:
        F1, F2, F3 = 32, 16, 28
    elif len(ft) == 2:
        F1, F2, F3 = ft[0], ft[1], 10
    elif len(ft) == 3:
        F1, F2, F3 = ft[0], ft[1], ft[2]
    # 4.input values: (x, y, z) 
    x_mesh, y_mesh = np.meshgrid(x, y)

    # make the meshgrid more intense. 

    # 5. Select the contour plot color.
    if MyColor == "PiYG": # @ "PiYG" color.
        contour = plt.contourf(x_mesh, y_mesh, z, cmap="PiYG", levels=400)
        custom_cmap = MyColor
    elif MyColor == "Accent":
        contour = plt.contourf(x_mesh, y_mesh, z, cmap="Accent", levels=400)
        custom_cmap = MyColor
    elif MyColor == "set3":
        contour = plt.contourf(x_mesh, y_mesh, z, cmap="set3", levels=400)
        custom_cmap = MyColor
    elif MyColor == "YlOrBr": # @ "YlOrBr" color.
        contour = plt.contourf(x_mesh, y_mesh, z, cmap="YlOrBr", levels=20)
        custom_cmap = MyColor
    elif MyColor == "jet":
        contour = plt.contourf(x_mesh, y_mesh, z, cmap="jet", levels=200)
        custom_cmap = MyColor
    elif MyColor == "viridis":
        contour = plt.contourf(x_mesh, y_mesh, z, cmap="viridis", levels=200)
        custom_cmap = MyColor

    elif MyColor == 1: # @ France flag color.
        # colomap_1. Define the colormap according to the energy limitation. 
        # @According to the energy limitation.
        color_1 = [[1, 0.3, 0.3], [0.95, 0.95, 0.95]]
        color_2 = [[0.95, 0.95, 0.95], [0.3, 0.8, 1]]
        """Here is the color with pink to blue. """
        # color_1 = [[1, 0, 0], [0.95, 0.95, 0.95]] # Pink color. 
        
        if clim is None:
            _color_1 = LinearSegmentedColormap.from_list('color1', color_1, N=40)
            _color_2 = LinearSegmentedColormap.from_list('color2', color_2, N=40)
        else:
            _min = np.abs(clim[0]) # @According to the energy limitation.
            _max = np.abs(clim[1])
            
            _color_1 = LinearSegmentedColormap.from_list('color1', color_1, N=_min)
            _color_2 = LinearSegmentedColormap.from_list('color2', color_2, N=_max)
        list_color1 = [_color_1(i) for i in range(_color_1.N)] 
        list_color2 = [_color_2(i) for i in range(_color_2.N)]
        list_color = list_color1 + list_color2
        custom_cmap = LinearSegmentedColormap.from_list('FranceFlag', list_color)
        # colormap_2 Plot the contour figure.
        if clim is None:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=len(list_color))
        else:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=len(list_color), vmin=clim[0], vmax=clim[1]) 
    elif MyColor == "1MWB": # Pink whit blue
        # colomap_1. Define the colormap according to the energy limitation. 
        # @According to the energy limitation.
        color_1 = [[0.8, 0, 0.8], [0.95, 0.95, 0.95]]
        color_2 = [[0.95, 0.95, 0.95], [0.3, 0.8, 1]]
        """Here is the color with pink to blue. """
        # color_1 = [[1, 0, 0], [0.95, 0.95, 0.95]] # Pink color. 
        
        if clim is None:
            _color_1 = LinearSegmentedColormap.from_list('color1', color_1, N=40)
            _color_2 = LinearSegmentedColormap.from_list('color2', color_2, N=40)
        else:
            _min = np.abs(clim[0]) # @According to the energy limitation.
            _max = np.abs(clim[1])
            
            _color_1 = LinearSegmentedColormap.from_list('color1', color_1, N=_min)
            _color_2 = LinearSegmentedColormap.from_list('color2', color_2, N=_max)
        list_color1 = [_color_1(i) for i in range(_color_1.N)] 
        list_color2 = [_color_2(i) for i in range(_color_2.N)]
        list_color = list_color1 + list_color2
        custom_cmap = LinearSegmentedColormap.from_list('FranceFlag', list_color)
        # colormap_2 Plot the contour figure.
        if clim is None:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=len(list_color))
        else:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=len(list_color), vmin=clim[0], vmax=clim[1]) 

    elif MyColor == 2:  
        # The map change for in a specific change.  
        color_1 = [[1, 0, 0], [0.95, 0.95, 0.95]]  
        color_2 = [[0.95, 0.95, 0.95], [0, 151/256, 244/256]]  
        color_3 = [[0, 151/256, 244/256], [0, 0, 1]]  
        if clim is None: 
            _min = 10  # @According to the energy limitation.  
            _middle = 20  
            _max = 170  
        else:
            _unit = (clim[1] - clim[0]) / 200
            _min = int((clim[2]- clim[0])/_unit) 
            _middle = int((clim[3]- clim[2])/_unit) 
            _max = int((clim[4]- clim[3])/_unit)
        _color_1 = LinearSegmentedColormap.from_list('color1', color_1, N=_min)  
        _color_2 = LinearSegmentedColormap.from_list('color2', color_2, N=_middle)  
        _color_3 = LinearSegmentedColormap.from_list('color2', color_3, N=_max)  
        list_color1 = [_color_1(i) for i in range(_color_1.N)]  
        list_color2 = [_color_2(i) for i in range(_color_2.N)]  
        list_color3 = [_color_3(i) for i in range(_color_3.N)]  
        list_color = list_color1 + list_color2 + list_color3  
        custom_cmap = LinearSegmentedColormap.from_list('FranceFlag', 
                                                        list_color)  
        # colormap_2 Plot the contour figure.  
        if clim is None:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, 
                                   levels=len(list_color))
        else:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, 
                                   levels=len(list_color), vmin=clim[0], 
                                   vmax=clim[1])
    elif MyColor == "Gradient":
        color = np.array([[16, 70, 128], [49, 124, 183], [109, 173, 209], 
                          [182, 215, 232],[233, 241, 244]])/256 
        
        color = np.flipud(color)
        # Smooth the color map. 
        custom_cmap = LinearSegmentedColormap.from_list("smooth_cmap", color, 
                                                        N=24) 
        if clim is None:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=24)

        else:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, 
                                   vmin=clim[0], vmax=clim[1], levels=70)
              
    else:
        # colormap_3 Define the colormap with average change.
        Colors = [(1, 0, 0), (1, 1, 0.9), (0.392, 0.584, 0.929)] # Cornflower blue to white to red
        cmap_name = 'custom_gradient'
        custom_cmap = LinearSegmentedColormap.from_list(cmap_name, Colors, N=100)
        if clim is None:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=100)
        else:
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, levels=100, vmin=clim[0], vmax=clim[1])

    # 5. control the contour line _con_line_.  
    if con_line is None:
        pass 
    else:
        if con_line[2] == 0:
            if clim is None:
                contour_lines = plt.contour(x_mesh, y_mesh, z, colors='black', linewidths=0.5)
            else: 
                contour_levels1 = list(np.linspace(clim[0], 0, 6))
                contour_levels2 = list(np.linspace(0, clim[1], 6))
                contour_levels = contour_levels1 + contour_levels2[1:]
                contour_lines = plt.contour(x_mesh, y_mesh, z, colors='black', levels=contour_levels, linewidths=0.5, vmin=clim[0], vmax=clim[1]) # @ contour lines
        else:
            if clim is None:
                clim = con_line 
            con_levels = np.linspace(con_line[0], con_line[1], con_line[2])
            contour_lines = plt.contour(x_mesh, y_mesh, z, levels=con_levels, 
                                        colors='black', linewidths=0.5) # @ contour lines
            contour = plt.contourf(x_mesh, y_mesh, z, cmap=custom_cmap, 
                                   vmin=clim[0], vmax=clim[1], levels=con_levels)

        def custom_format(val):
            if ContourLable ==None:
                return f'${val:.1f}$'
            elif ContourLable =="10^-val":
                return f'$10^{{{-val:.1f}}}$'
            else:
                return f'${val:.1f}$'

        plt.clabel(contour_lines, inline=True, fontsize=14, fmt=custom_format)

    # 6. Set x and y axis scale
    if xylog is None:  
        pass  
    elif xylog in ["xlog", "Xlog"]:  
        plt.xscale('log')  
    elif xylog in ["ylog", "Ylog"]:  
        plt.yscale('log')  
    else:  
        plt.xscale('log')  
        plt.yscale('log')

    plt.xlabel(xname, fontsize=F1, fontname="Arial", weight='bold')  
    plt.ylabel(yname, fontsize=F1, fontname="Arial", weight='bold')  
    plt.title(title, fontsize=F1, fontname="Arial", weight='bold')

    #------------------------------------------------------------------------
    # 7. Set the label, ticks, and grid.
    if xticks is None:  
        plt.xticks(weight='bold')  
    else:  
        x_tick = xticks["xticks"]  
        x_labels = xticks["xlables"]  
        plt.xticks(x_tick, x_labels, weight='bold')  
    if yticks is None:  
        plt.yticks(weight='bold')  
    else:  
        y_tick = yticks["yticks"]  
        y_labels = yticks["ylables"]  
        plt.yticks(y_tick, y_labels,  weight='bold')

    plt.tick_params(axis='both',       # apply to both x and y axes
                which='major',     # apply to major ticks
                direction='out',    # 'in', 'out', or 'inout'
                length=6,          # length of the ticks
                width=2,           # width of the ticks
                labelsize=F2)      # font size of tick labels
    
    if grid is not None:  
        plt.grid(True, which='both', axis='both', linestyle='--', color="black")
    
    #------------------------------------------------------------------------
    # 8. Color bar. 
    if pltcbar == 1:  
        cbar = plt.colorbar(contour) # Create color bar
        # Adjust font size and weight for color bar label
        cbar.set_label(cbar_name, fontsize=F3, fontweight='bold') 
        cbar.ax.tick_params(labelsize=F2) # Adjust font size for color bar ticks
        if cticks is None:  
            cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%1.2f'))  
        else:  
            cbar_ticks = cticks["cticks"]  
            ctick_labels = cticks["clables"]  
            cbar.set_ticks(cbar_ticks)  
            cbar.set_ticklabels(ctick_labels, fontweight="bold")

    #------------------------------------------------------------------------
    # 9. Plot lines in the contour plot: _ax_line_
    if ax_line is not None:
        plt.axvline(x=ax_line[0], color='black', linestyle='-', linewidth=2)
        plt.axhline(y=ax_line[1], color='black', linestyle='-', linewidth=2)

    #-----------------------------------------------------------------------
    if adjust is None:
        pass 
    else: 
        plt.subplots_adjust(top=adjust[0], left=adjust[1], bottom=adjust[2], 
                            right=adjust[3])
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)  # Change 2 to your desired width




# 0_2 Plot function: contour plot with size of (2,1). 
def subplot_cont(x1, y1, z1, label1, x2, y2, z2, label2, subfig=None, MyColor=0, xylog=None, ft=None, clim1=None, clim2=None, con_lim1=None, ax_line=None):
    if subfig is None:
        plt.figure(figsize=(16, 10))
    else:
        plt.figure(figsize=(subfig[0], subfig[1]))
    # Plot 1: E_binding
    plt.subplot(1, 2, 1)
    contour_plot(x1,y1,z1,label=label1, pltfig=[1], MyColor=MyColor, xylog=xylog, ft=ft, clim=clim1, con_lim=con_lim1, ax_line=ax_line)
    plt.tight_layout(pad=2.0, w_pad=8, h_pad=1.0)
    # Plot 2: theta
    plt.subplot(1, 2, 2)
    contour_plot(x2,y2,z2,label=label2, pltfig=[1], MyColor=MyColor, xylog=xylog, ft=ft, clim=clim2, con_lim=con_lim1, ax_line=ax_line)


# 0_3 Plot function: contour plot with size of (2,2). 
def subplot_cont22(x1, y1, zs, labels, subfig=None, MyColor=0, xylog=None, ft=None, clims=None, con_lims=None, ax_line=None):
    if subfig is None:
        plt.figure(figsize=(16, 10))
    else:
        plt.figure(figsize=(subfig[0], subfig[1]))
    if con_lims is None:
        con_lims = [None, None, None, None]
    plt.subplot(2, 2, 1)
    contour_plot(x1,y1,zs[0],label=labels[0], pltfig=[1], MyColor=MyColor, xylog=xylog, ft=ft, clim=clims[0], con_lim=con_lims[0], ax_line=ax_line)
    plt.subplot(2, 2, 2)
    contour_plot(x1,y1,zs[1],label=labels[1], pltfig=[1], MyColor=MyColor, xylog=xylog, ft=ft, clim=clims[1], con_lim=con_lims[1], ax_line=ax_line)
    plt.subplot(2, 2, 3)
    contour_plot(x1,y1,zs[2],label=labels[2], pltfig=[1], MyColor=MyColor, xylog=xylog, ft=ft, clim=clims[2], con_lim=con_lims[2], ax_line=ax_line)
    plt.subplot(2, 2, 4)
    contour_plot(x1,y1,zs[3],label=labels[3], pltfig=[1], MyColor=MyColor, xylog=xylog, ft=ft, clim=clims[3], con_lim=con_lims[3], ax_line=ax_line)
    plt.tight_layout(pad=2.0, w_pad=8, h_pad=1.0)