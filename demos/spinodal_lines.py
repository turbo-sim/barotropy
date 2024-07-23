    # # Create figure for plotting
    # if i == 0:
    #     # Plot phase diagram
    #     fig_1, ax_1 = plt.subplots(figsize=(6.0, 5.0))
    #     ax_1.set_xlabel("Entropy (J/kg/K)")
    #     ax_1.set_ylabel("Temperature (K)")
    #     # s_min, s_max = 1.3, 1.7
    #     # T_min, T_max = 25, 35
    #     # ax.set_xlim([s_min, s_max])
    #     # ax.set_ylim([T_min, T_max])
    #     ax_1 = fluid.plot_phase_diagram(
    #         var_x,
    #         var_y,
    #         axes=ax_1,
    #         plot_critical_point=True,
    #         plot_saturation_line=True,
    #         plot_spinodal_line=False,
    #         plot_quality_isolines=True,
    #     )

    #     # Plot density vs pressure
    #     fig_2, ax_2 = plt.subplots(figsize=(6.0, 5.0))
    #     ax_2.set_xlabel("Pressure (Pa)")
    #     ax_2.set_ylabel("Density (kg/m$^3$)")
