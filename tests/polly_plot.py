from pybench import parser, Benchmark
import pylab as plt


class PollyPlot(Benchmark):
    """
    example: python tests/polly_plot.py -ts "2 4 8" "4 4 1000" "default" "nopolly"
    """

    figsize = (6, 4)


if __name__ == '__main__':
    p = parser(description="polly performance benchmark for different tile sizes ")
    p.add_argument('-ts', '--tile_size', type=str, nargs='+',
                   help='tile sizes to plot')

    args = p.parse_args()
    ts = args.tile_size
    b = PollyPlot(benchmark='polly_tile_bench', resultsdir=args.resultsdir, plotdir=args.plotdir)

    b.combine_series([('tile_size', ts)], filename='polly_tile_bench')

    res = []
    for data in b.result['timings'].values():
        del(data['total'])
        res.append(data.values()[0])

    labels_ = b.result['timings'].keys()
    labels = [x[0] for x in labels_]

    fig = plt.figure()
    ax = fig.add_axes((.1, .4, .8, .5))

    plt.bar(range(0, len(res)), res, width=0.5, align='center')
    plt.xticks(range(0, len(res)), labels, rotation=17)

    plt.ylabel('running time ')
    plt.savefig('tile_test_plot')
