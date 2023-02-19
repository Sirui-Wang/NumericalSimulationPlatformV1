import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
from tqdm import tqdm


def update1(frame_key):
    global SaveDict, LineLevy, LineNormal, x_text, LineHist, LineNormPDF, LineLevyPDF, loc, scale, ax1, ax2
    x = SaveDict[frame_key]["x"]
    x_norm= SaveDict[frame_key]["x_norm"]
    ASSeries= SaveDict[frame_key]["ASSeries"]
    NSeries= SaveDict[frame_key]["NSeries"]
    Normal_PDF= SaveDict[frame_key]["Normal_PDF"]
    Levy_PDF= SaveDict[frame_key]["Levy_PDF"]
    LineLevy.set_data(np.arange(0,len(ASSeries)),ASSeries)
    LineNormal.set_data(np.arange(0,len(NSeries)),NSeries)
    LineNormPDF.set_data(x_norm, Normal_PDF)
    LineLevyPDF.set_data(x, Levy_PDF)
    x_text.set_text("alpha, beta = {}".format(frame_key))
    ax1.set_xlim([0,len(ASSeries)])
    ax2.set_xlim([x_norm[0],x_norm[-1]])
    ax1.set_ylim([min(ASSeries),max(ASSeries)])
    ax2.set_ylim([min(Levy_PDF), max(Normal_PDF)])
    return LineLevy, LineNormal, LineNormPDF, LineLevyPDF, ax1, ax2



from matplotlib import animation
from matplotlib.animation import PillowWriter


def main():
    ss.levy_stable.parameterization = "S1"
    global SaveDict, LineLevy, LineNormal, x_text, LineHist, LineNormPDF, LineLevyPDF, loc, scale,ax1, ax2
    loc = 2 # similar to mean in normal distribution
    scale = 1 # similar to std in normal distribution
    change_alpha = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
    # change_alpha = [1.0,1.5]
    change_beta = [-0.5]
    SaveDict = {}
    boundary = 10
    change_step = 0.1
    x = np.arange(-boundary,boundary+change_step,change_step)
    for i in tqdm(range(len(change_beta))):
        for j in tqdm(range(len(change_alpha))):
            alpha = change_alpha[j]
            beta = change_beta[i]
            ASSeries = ss.levy_stable.rvs(alpha, beta, loc=loc, scale=scale*1/np.sqrt(2), size = 100)
            NSeries = ss.norm.rvs(loc=loc, scale=scale, size=1000)
            Normal_PDF = ss.norm.pdf(x, loc=loc, scale=scale)
            Levy_PDF = ss.levy_stable.pdf(x, alpha, beta, loc=loc, scale=scale*1/np.sqrt(2))
            SaveDict[(alpha, beta)] = {"x":x,"x_norm":x, "ASSeries":ASSeries,"NSeries":NSeries,"Normal_PDF":Normal_PDF,"Levy_PDF":Levy_PDF}
    fig, (ax1, ax2) = plt.subplots(2,1)
    fig.suptitle("Animated Levy Distribution")
    LineLevy, = ax1.plot([],"r-", label="levy_stable noise time series")
    LineNormal, = ax1.plot([],"b-", label="normal noise time series")
    x_text = ax2.text(-3,0,"", fontsize = 15)
    # LineHist, = ax2.hist([], density=True, bins='auto', histtype='stepfilled', alpha=0.8)
    LineNormPDF, = ax2.plot([], [],'b-', lw=5, alpha=0.6, label='normal pdf')
    LineLevyPDF, = ax2.plot([], [],'r-', lw=1, alpha=1, label='levy_stable pdf')
    ax1.legend(loc='best', frameon=False)
    ax2.legend(loc='best', frameon=False)
    Figureanimation1 = animation.FuncAnimation(fig, update1, interval=500, frames=list(SaveDict.keys()), blit=False)
    Figureanimation1.save("Levy_Distribution.gif", dpi=300, writer=PillowWriter(fps=1))
    plt.show()


main()


# alpha = 1  # characteristic exponents that discribes the tail of the distribution. （0，2] -- 0<alpha<=2
# beta = -1 # skewness that defines if the distribution is left skewed (beta <0) or right skewed (beta >0) [-1,1]
# loc = 0 # similar to mean in normal distribution
# scale = 1 # similar to std in normal distribution
# # x = np.linspace(ss.levy_stable.ppf(-0.99, alpha, beta,loc=loc, scale=scale),
# #                 ss.levy_stable.ppf(0.99, alpha, beta,loc=loc, scale=scale), 100)
# # x_norm = np.linspace(ss.norm.ppf(0.001, loc=loc, scale=scale*np.sqrt(2)), ss.norm.ppf(0.999,loc=loc, scale=scale*np.sqrt(2)),100)
# x = np.arange(-4,4.01,0.1)
# ASSeries = ss.levy_stable.rvs(alpha, beta, loc=loc, scale=scale, size = 1000)
# NSeries = ss.norm.rvs(loc=loc, scale=scale*np.sqrt(2), size=1000)
# fig, ax = plt.subplots(2,1)
# ax[0].plot(ASSeries,"r-", label="levy_stable noise time series")
# ax[0].plot(NSeries,"b-", label="normal noise time series")
# ax[1].hist(ASSeries, density=True, bins="auto", histtype='stepfilled', alpha=0.5)
# ax[1].set_xlim([x[0], x[-1]])
# ax[1].plot(x, ss.norm.pdf(x, loc=loc, scale=scale*np.sqrt(2)),
#         'b-', lw=5, alpha=0.6, label='normal pdf')
# ax[1].plot(x, ss.levy_stable.pdf(x, alpha, beta, loc=loc, scale=scale),
#         'r-', lw=1, alpha=1, label='levy_stable pdf')
# ax[0].legend(loc='best', frameon=False)
# ax[1].legend(loc='best', frameon=False)
# plt.show()
