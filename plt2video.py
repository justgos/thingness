import matplotlib.pyplot as plt
import numpy as np


def fig2data ( fig ):
    """
    @brief Convert a Matplotlib figure to a 4D numpy array with RGBA channels and return it
    @param fig a matplotlib figure
    @return a numpy 3D array of RGBA values
    Ref: http://www.icare.univ-lille1.fr/tutorials/convert_a_matplotlib_figure
    """
    # Get rid of the margins
    fig.gca().set_axis_off()
    fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
                hspace = 0, wspace = 0)
    #fig.margins(0,0)
    from matplotlib.ticker import NullLocator
    fig.gca().xaxis.set_major_locator(NullLocator())
    fig.gca().yaxis.set_major_locator(NullLocator())
    
    # draw the renderer
    fig.canvas.draw ( )

    # Get the RGBA buffer from the figure
    w,h = fig.canvas.get_width_height()
    abuf = np.frombuffer(fig.canvas.tostring_argb(), dtype=np.uint8).reshape((h, w, 4))[:, :, :1]
    buf = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8).reshape((h, w, 3))
    buf = (buf / 255.0 * abuf + ((255-abuf) / 255.0) * np.array([[[255, 255, 255]]])).astype(np.uint8)
    #buf = buf.transpose((1, 0, 2))
    #buf = np.fromstring ( fig.canvas.tostring_argb(), dtype=np.uint8 )
    #buf.shape = ( w, h,4 )
    
    # canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
    #buf = np.roll ( buf, 3, axis = 2 )
    return buf


def frames2web(frames, fps=60):
    """Converts NxHxWxC NumPu array into Base64-encode video, wrapped in HTML (for Visdom visualization)"""
    import skvideo.io
    import tempfile
    import mimetypes
    import base64
    videofile = '/tmp/%s.mp4' % next(tempfile._get_candidate_names())
    skvideo.io.vwrite(videofile, frames,
                      inputdict={
                        "-r": "%d" % fps
                      },
                      outputdict={
                        "-b": "4000000"
                      })

    extension = videofile[-3:].lower()
    mimetypes = dict(mp4='mp4', ogv='ogg', avi='avi', webm='webm')
    mimetype = mimetypes.get(extension)

    with open(videofile, 'rb') as f:
        bytestr = f.read()
    videodata = """
        <video width="%d" controls>
            <source type="video/%s" src="data:video/%s;base64,%s">
            Your browser does not support the video tag.
        </video>
    """ % (frames.shape[2], mimetype, mimetype, base64.b64encode(bytestr).decode('utf-8'))
    return videodata

from IPython.display import HTML


from visdom import Visdom
viz = Visdom()
viz_opts = {
    'sim_slice': dict()
}
viz_windows = { k: None for k in viz_opts.keys() }
# vizdom .video() modified to use scikit-video instead of OpenCV
def viz_video(frames, win=None, env=None, opts=None, fps=60):
    return viz.text(text=frames2web(frames, fps), win=win, env=env, opts=opts)

