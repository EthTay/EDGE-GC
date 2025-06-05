import matplotlib.pyplot as plt
import pynbody
import numpy as np 
import pandas as pd
import gc

import json

from plot_data import load_data
from matplotlib.colors import LinearSegmentedColormap

def load_group_data(sim):
    path = "/user/HS501/et00737/starcluster_Stuff/%s.starcluster_data" % sim
    loaded_data = [[], []]
    with open(path, "r") as f:
        total_groups = int(f.readline())
        for i in range(total_groups):
            loaded_data[0].append(int(f.readline()))
            num_part = int(f.readline())
            particle_data = np.zeros((num_part, 2), dtype=int)
            for j in range(num_part):
                read_buffer = f.readline().split()
                particle_data[j, 0] = int(read_buffer[0])
                particle_data[j, 1] = int(read_buffer[1])
            loaded_data[1].append(particle_data)
    return loaded_data

def get_snap_info(sim,internal_index):

  snap_data = load_group_data(sim)

  snap = snap_data[0][internal_index]
  star_IDs = snap_data[1][internal_index][snap_data[1][internal_index][:, 1] == 4][:, 0]
  DM_IDs = snap_data[1][internal_index][snap_data[1][internal_index][:, 1] == 1][:, 0]

  return sim,snap,star_IDs,DM_IDs

def mutual(halo_ids,check_ids):
   
  tf=np.in1d(halo_ids,check_ids)

  return (len(tf[tf])**2)/(len(halo_ids)*len(check_ids))

def main(gas_subs,dm_subs):
  with open("/vol/ph/astro_data/shared/etaylor/cmaps/EDGE_cmap.json", "r") as f:
      EDGE_cmap = LinearSegmentedColormap("/vol/ph/astro_data/shared/etaylor/cmaps/EDGE_cmap.json", json.load(f))
  
  GC_df = pd.read_csv("/user/HS501/et00737/starcluster_Stuff/data/refresh.csv")
  dwarfs_df = pd.read_csv('/user/HS501/et00737/starcluster_Stuff/data/refresh_dwarfs.csv')
  
  observational_dwarfs = load_data(1)
  observational_GC = load_data(0)
  
  target_sims = ['Halo1445_fiducial_hires']#,'Halo605_fiducial_hires','Halo624_fiducial_hires','Halo1445_fiducial_hires','Halo1459_fiducial_hires']
  target_index = [2]#,2,1,2,0] 
  titles = ['DMC mini-halo Formation']#,'Gas Precipitation\nvia merger','Triggered Formation\nvia merger','DMC mini-halo formation','DMC mini-halo Formation']
  
  EDGE_path = '/vol/ph/astro_data/shared/morkney/EDGE/'
  
  for i,buffer_data in enumerate(zip(target_sims,target_index)):
    simulation,index = buffer_data
    mask = ((GC_df['Simulation']==simulation) &(GC_df['Internal ID']==index))
    snap=list(GC_df['Output Number'][mask])[0]
  
    half_light_radius = list(GC_df['Half-light radius'][mask])[0]
  
  
    buffer_sim,buffer_snap,star_IDs,DM_IDs = get_snap_info(simulation,index)
  
    s = pynbody.load('%s%s/output_%s' %(EDGE_path,simulation,f'{snap:05d}'))
    s.physical_units('pc')
    stars_tf = np.in1d(s.s['iord'],star_IDs)
    h=s.halos()
    halo_star_IDs = h[1].s['iord']
    halo_dm_IDs = h[1].d['iord']
  
    cen=pynbody.analysis.halo.shrink_sphere_center(h[1])
    s['pos']-=cen
    x_mean,y_mean,z_mean=s.s[stars_tf]['x'].mean(),s.s[stars_tf]['y'].mean(),s.s[stars_tf]['z'].mean()
  
    time=[]
    cen_halo_ids=[]
  
    for j in range(-2,3):
       circle=plt.Circle((x_mean, y_mean), half_light_radius, color="w", fill=False,ls='--')
       s=None
       h=None
       gc.collect()
       s = pynbody.load('%s%s/output_%s' %(EDGE_path,simulation,f'{snap+j:05d}'))
       s.physical_units('pc')
       h=s.halos()
       if j <0:
         mutual_c = np.array([mutual(h[h_index+1].d['iord'],DM_IDs) for h_index in range(len(h))])
         host_index = mutual_c.argsort()[-1]+1
         print(mutual_c[mutual_c.argsort()[-1]])
         cen=pynbody.analysis.halo.shrink_sphere_center(h[host_index])
       elif j>=0:
         stars_tf = np.in1d(s.s['iord'],star_IDs)
         x_mean,y_mean,z_mean=s.s[stars_tf]['x'].mean(),s.s[stars_tf]['y'].mean(),s.s[stars_tf]['z'].mean()
         cen=[x_mean,y_mean,z_mean]
      
       s['pos']-=cen
       stars_tf = np.in1d(s.s['iord'],star_IDs)
       time.append(s.properties['time'].in_units('Gyr'))
       cen_halo_ids.append(host_index)
    
  
  ###plot another
  
       pynbody.plot.image(s.g[s.g['r'] <= 6000], units='Msol kpc^-2',cmap='viridis',width = '4000 pc' ,show_cbar=False,subplot=gas_subs[2+j]) #1459 0 
       gas_subs[2+j].scatter(s.s["x"][~stars_tf], s.s["y"][~stars_tf], c="black", s=5, alpha=0.1)
       gas_subs[2+j].scatter(s.s[stars_tf]["x"], s.s[stars_tf]["y"], c="red", s=5, alpha=0.2)
    
  ###plot the last
  
       pynbody.plot.image(s.d[s.d['r'] <= 6000], units='Msol kpc^-2',cmap=EDGE_cmap,width = '4000 pc' ,show_cbar=False,subplot=dm_subs[2+j]) #1459 0
       dm_subs[2+j].scatter(s.s["x"][~stars_tf], s.s["y"][~stars_tf], c="black", s=5, alpha=0.1)
       dm_subs[2+j].scatter(s.s[stars_tf]["x"], s.s[stars_tf]["y"], c="red", s=5, alpha=0.2)
    
  
       gas_subs[j+2].set_ylim(-1750, 1750)
       gas_subs[j+2].set_xlim(-1750, 1750)
       dm_subs[j+2].set_ylim(-1750, 1750)
       dm_subs[j+2].set_xlim(-1750, 1750)
  
  
  time=np.array(time)*1000
  time-=list(GC_df["Mean Tform"][mask])[0]
  
  for i in range(len(time)):
      gas_subs[i].set_xlabel("")
      dm_subs[i].set_title("")
  
  for i in range(1,len(time)):
    gas_subs[i].set_ylabel("")
    dm_subs[i].set_ylabel("")
  
  plt.ion()
  plt.show()
  
  target_sims = ['Halo1445_fiducial_hires']#,'Halo605_fiducial_hires','Halo624_fiducial_hires','Halo1445_fiducial_hires','Halo1459_fiducial_hires']
  target_index = [2]#,2,1,2,0]
  titles = ['DMC mini-halo Formation']#,'Gas Precipitation\nvia merger','Triggered Formation\nvia merger','DMC mini-halo formation','DMC mini-halo Formation']
  
  position = 725/750
  left = -700/750
  
  for j in range(len(gas_subs)):
      for i in range(len(gas_subs[j].texts)):
          gas_subs[j].texts[i].set_visible(False)
  
  
  gas_subs[0].text(
      abs(left)*gas_subs[0].get_xlim()[0],
      position*gas_subs[0].get_ylim()[1],
      "%s" % (titles[0]),
      color="white",
      horizontalalignment="left",
      verticalalignment="top",
      size=10,
  )
  
  titles = ['DMC mini-halo Formation']
  
  for i in range(2):#len(time)):
      gas_subs[i].text(position*gas_subs[0].get_xlim()[1], position*gas_subs[0].get_ylim()[1], "T\N{MINUS SIGN}%s Myr" % (abs(np.round(time[i], 1))), color="white",horizontalalignment='right',verticalalignment="top",weight='bold',size=13)
  
  for i in range(2,len(time)):
      #gas_subs[i].text(position*gas_subs[0].get_xlim()[0], position*gas_subs[0].get_ylim()[1], "T\N{PLUS SIGN}%s Myr" % (abs(np.round(time[i], 1))), color="white",horizontalalignment='right',verticalalignment="top",weight='bold',size=13)
      gas_subs[i].text(position*gas_subs[0].get_xlim()[1], position*gas_subs[0].get_ylim()[1], "T\N{PLUS SIGN}%s Myr" % (abs(np.round(time[i], 1))), color="white",horizontalalignment='right',verticalalignment="top",weight='bold',size=13)
  
  
  dm_subs[0].set_yticks(
      [-1500,-1000, -500, 0,500,1000,1500],
      [
          "\N{MINUS SIGN}1500",
          "\N{MINUS SIGN}1000",
          "\N{MINUS SIGN}500",
          "0",
          "500",
          "1000",
          "1500"
      ],size=13
  )
    
  gas_subs[0].set_yticks(
      [-1500,-1000, -500, 0,500,1000,1500],
      [
          "\N{MINUS SIGN}1500",
          "\N{MINUS SIGN}1000",
          "\N{MINUS SIGN}500",
          "0",
          "500",
          "1000",
          "1500"
      ],size=13
  )
  
  for i in range(len(time)):
    
    dm_subs[i].set_xticks(
        [],
        [      ],size=13
    )
    
    dm_subs[i].set_xlabel(dm_subs[i].get_xlabel(),size=20)
  
  for i in range(len(time)):
    gas_subs[i].tick_params("x", direction="in")
  for i in range(1,len(time)):
    gas_subs[i].tick_params("y", direction="in")
    dm_subs[i].tick_params("y", direction="in")
  
  gas_subs[0].set_ylabel("$y/\\mathrm{pc}$", size=17.5)
  dm_subs[0].set_ylabel("$y/\\mathrm{pc}$", size=17.5)
  
#  for i in range(20):
#    gs.tight_layout(plt.gcf(), h_pad=0, w_pad=0)
  
  for i in range(len(time)):
  
      dm_subs[i].set_xticks(
          [-1500,-1000,-500, 0,500,1000,1500],
          [   "",
              "\N{MINUS SIGN}1000",
              "\N{MINUS SIGN}500",
              "0",
              "500",
              "1000",
              ""
          ],size=13
          )
  
      dm_subs[i].set_xlabel(dm_subs[i].get_xlabel(), size=20)
  
  
#  gs.update()
  
  
  
  
  
  for j in range(len(gas_subs)):
      for i in range(len(gas_subs[j].texts)):
          gas_subs[j].texts[i].set_visible(False)
      for i in range(len(gas_subs[j].lines)):
          gas_subs[j].lines[i].set_visible(False)
  text_left = gas_subs[0].get_xlim()[0]*(725/750)
  text_top = gas_subs[0].get_ylim()[1]*(725/750)
  text_bottom = gas_subs[0].get_ylim()[0]*(725/750)
  text_right = gas_subs[0].get_xlim()[1]*(725/750)
  gas_subs[1].text(
      -(700/750)*1750,
      -(725/750)*1750,
      "%s" % (titles[0]),
      color="white",
      horizontalalignment="left",
      verticalalignment="bottom",
      size=25,
  )
  for i in range(2):
      gas_subs[i].text(text_right, text_top, "T\N{MINUS SIGN}%s Myr" % (abs(np.round(time[i], 1))), color="white",horizontalalignment='right',verticalalignment="top",weight='bold',size=20)
  for i in range(2,len(time)):
      gas_subs[i].text(text_right, text_top, "T\N{PLUS SIGN}%s Myr" % (abs(np.round(time[i], 1))), color="white",horizontalalignment='right',verticalalignment="top",weight='bold',size=20)
  for i in range(len(gas_subs)):
    dm_subs[i].set_xticks(
         [],
         [],size=13
     )
  dm_subs[0].set_yticks([],[])
  gas_subs[0].set_yticks([],[])
  for i in range(len(time)):
    dm_subs[i].set_xticks([],[])
    dm_subs[i].set_xlabel("")
    dm_subs[i].set_xticks(
         [-600, -300,  0, 300, 600],
         [
             "\N{MINUS SIGN}600",
             "\N{MINUS SIGN}300",
             "0",
             "300",
             "600",
         ],size=13
     )
    dm_subs[i].set_xlabel(dm_subs[i].get_xlabel(),size=20)
  for i in range(len(time)):
    gas_subs[i].tick_params("x", direction="in")
  for i in range(1,len(time)):
    gas_subs[i].tick_params("y", direction="in")
    dm_subs[i].tick_params("y", direction="in")
  gas_subs[0].set_ylabel("$y/\\mathrm{pc}$", size=17.5)
  dm_subs[0].set_ylabel("$y/\\mathrm{pc}$", size=17.5)
  gas_subs[0].set_ylabel("$y/\\mathrm{pc}$", size=17.5)
  dm_subs[0].set_ylabel("$y/\\mathrm{pc}$", size=17.5)
  gas_subs[0].set_ylabel("", size=17.5)
  dm_subs[0].set_ylabel("", size=17.5)
  bar_left_position = gas_subs[1].get_xlim()[0]*(675/750)
  bar_right_position=bar_left_position+500
  text_left_position = gas_subs[0].get_xlim()[0]*(700/750)
  text_height_pos = (645/750)*gas_subs[0].get_ylim()[1]
  gas_subs[1].plot([bar_left_position, bar_right_position], [1400, 1400], "w", linewidth=4)
  gas_subs[1].text(
      text_left_position,
      text_height_pos,
      "500 pc",
      color="white",
      weight="bold",
      size=25,
      horizontalalignment="left",
  )
  
  for i in range(-1,1):
    gas_subs[i].set_visible(False)
    dm_subs[i].set_visible(False)
  
  return 0

#plt.savefig('/user/HS501/et00737/starcluster_Stuff/paper1/pdf/%s_%s_tight.pdf' %(target_sims[0],target_index[0]), bbox_inches = 'tight',pad_inches = 0)
#plt.savefig('/user/HS501/et00737/starcluster_Stuff/paper1/png/%s_%s_tight.png' %(target_sims[0],target_index[0]), bbox_inches = 'tight',pad_inches = 0)



