require 'Tioga/FigureMaker'
require './plot_styles.rb'
require 'Dobjects/Function'

class MyPlots

  include Math
  include Tioga
  include FigureConstants
  include MyPlotStyles

  def t
    @figure_maker
  end

  def initialize
    @figure_maker = FigureMaker.default

    t.tex_preview_preamble += "\n\t\\usepackage{mathtools}\n"
    t.tex_preview_preamble += "\n\t\\usepackage{amsmath}\n"
    t.tex_preview_preamble += "\n\t\\usepackage[charter]{mathdesign}\n"

    t.save_dir = 'plots'

    t.def_figure('conduction-convergence-study') do
      mnras_style
      enter_page_wide

      convergence_plot
    end

    t.def_figure('cond-test-temp-profiles') do
      mnras_style
      enter_page_single

      profile_plot
    end

    t.def_figure('cond-test-convergence') do
      mnras_style
      enter_page_single

      kappa_plot
    end

  end

  # square plot
  #
  def enter_page_single
    mnras_style

    t.xlabel_shift = 2.0
    t.ylabel_shift = 1.75

    t.default_frame_left   = 0.15
    t.default_frame_right  = 0.96
    t.default_frame_top    = 0.96
    t.default_frame_bottom = 0.15

    t.default_page_width  = 72 * 3.0

    t.default_page_height = t.default_page_width * \
      (t.default_frame_right - t.default_frame_left) / \
      (t.default_frame_top - t.default_frame_bottom)

    t.default_enter_page_function
  end

  # full page plot
  #
  def enter_page_wide
    mnras_style

    t.xlabel_shift = 2.0
    t.ylabel_shift = 1.75

    t.default_frame_left   = 0.07
    t.default_frame_right  = 0.98
    t.default_frame_top    = 0.96
    t.default_frame_bottom = 0.20

    t.default_page_width  = 72 * 7.0

    @margin = 0.2

    golden_ratio = 1.61803

    t.default_page_height = t.default_page_width * \
                            (t.default_frame_right - t.default_frame_left) / \
                            (t.default_frame_top - t.default_frame_bottom) / \
                            (2.0 + @margin) / golden_ratio

    t.default_enter_page_function
  end

  def convergence_plot
    space = (1.0 + @margin) / (2.0 + @margin)

    t.subfigure('right_margin' => space) { profile_plot }
    t.subfigure('left_margin'  => space) { kappa_plot }
  end


  def profile_plot
    t.do_box_labels(nil, '$r$', '$T$')

    t.yaxis_locations_for_major_ticks = [1.0, 1.025, 1.05, 1.075, 1.1]
    t.yaxis_tick_labels = %w{1.0 1.025 1.05 1.075 1.1}


    t.show_plot_with_legend('plot_right_margin' => 0.0,
                            'legend_left_margin' => 0.05) do
      t.show_plot([0.0, 1.0, 1.11, 1.0]) do
        r, t0, t1, t2 = Dvector.fancy_read('analytic/analytic.dat')

        t.show_polyline(r, t0, MidnightBlue,  '\Large $t = 0$')
        t.show_polyline(r, t1, FireBrick,     '\Large $t = t_{\text{cond}}/2$')
        t.show_polyline(r, t2, DarkGoldenrod, '\Large $t = t_{\text{cond}}$')

        old_width = t.line_width
        t.line_width = 0.75


        simdir = "simulation-results/numeric-profiles"

        _, x, y = Dvector.fancy_read("#{simdir}/256-end.dat")

        t.show_polyline(x, y, FireBrick, nil, Line_Type_Dash)

        _, x, y = Dvector.fancy_read("#{simdir}/128-end-long.dat")

        t.show_polyline(x, y, DarkGoldenrod, nil, Line_Type_Dash)

        t.line_width = old_width
      end
    end
  end


  def kappa_plot
    t.do_box_labels(nil, '$\sigma / \Delta x \sim (\Delta x \, \nabla \ln T)^{-1}$',
                   '$\kappa_{\perp} / \kappa_{\parallel}$')

    t.yaxis_log_values = true

    t.xaxis_locations_for_major_ticks = [0, 1, 2, 3, 4, 5, 6]
    t.xaxis_tick_labels = %w{1 2 4 8 16 32 64}

    t.show_plot([Math.log2(1), Math.log2(64), 1.log10, (1e-5).log10]) do

      # show the mathematica fit to the data
      #
      x = Dvector.new([0.0, Math.log2(64)])
      y = -1.44448 - 1.93562 * x

      old_width = t.line_width
      t.line_width = 0.75
      t.show_polyline(x/Math.log(2.0), y/Math.log(10.0),
                     MidnightBlue)
      t.line_width = old_width


      _, _, res2, kapparatio = Dvector.fancy_read('simulation-results/convergence.dat')

      kapparatio = -(kapparatio.safe_log10)
      res2 = (res2.safe_log10) / (2.log10)

      res2.each2(kapparatio) do |x,y|
        t.show_marker('x'      => x,
                      'y'      => y,
                      'marker' => Bullet,
                      'scale'  => 0.5)
      end
    end
  end

end

MyPlots.new

# Local Variables:
# fill-column: 96
# End:
