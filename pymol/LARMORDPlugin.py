""" 03_2016: Aaron T. Frank
      PyMOL plugin for running LARMORD and display chemical shift errors as spheres.

"""

# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
#               This PyMOL Plugin is Copyright (C) 2016 by 
#            Aaron T. Frank <afrankz at umich dot edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

# python lib
import os
import sys 
import platform
if sys.version_info >= (2,4):
    import subprocess # subprocess is introduced in python 2.4
import math
import random
import tempfile
import tkSimpleDialog
import tkMessageBox
import tkFileDialog
import tkColorChooser

import Tkinter

# pymol lib
try: 
    from pymol import cmd
    from pymol.cgo import *
    from pymol import stored
    from pymol import cmd
except ImportError:
    print 'Warning: pymol library cmd not found.'
    sys.exit(1)
    
# external lib    
try:
    import Pmw
except ImportError:
    print 'Warning: failed to import Pmw. Exit ...'
    sys.exit(1)
    
VERBOSE = False

#################
## here we go
#################
def __init__(self):
    """ LARMORD plugin for PyMol
    """
    self.menuBar.addmenuitem('Plugin', 'command','LARMORD', label = 'LARMORD', command = lambda s=self : LARMORDPlugin(s))
    

#################
## GUI related
#################
class LARMORDPlugin:
  
    def __init__(self, app):
        
        self.parent = app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons = ('Run LARMORD','Render Errors','Exit'),
                                 title = 'LARMORD Plugin for PyMOL',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        # parameters used by LARMORD
        self.pymol_sel     = Tkinter.StringVar()
        self.larmord_bin    = Tkinter.StringVar()
        self.larmord_para    = Tkinter.StringVar()
        self.larmord_ref    = Tkinter.StringVar()
        self.larmord_cs    = Tkinter.StringVar()
        self.larmord_rlt_dict = {}
        self.larmord_errors = {}
        self.larmord_error_color = Tkinter.StringVar()
        self.larmord_error_sel = Tkinter.StringVar()

        self.sel_obj_list = [] # there may be more than one seletion or object
                               # defined by self.pymol_sel
                               # treat each selection and object separately
        if 'LARMORD_BIN' not in os.environ and 'PYMOL_GIT_MOD' in os.environ:
            if sys.platform.startswith('linux') and platform.machine() == 'x86_32':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'],"Larmord","i86Linux2","larmord")
                os.environ['LARMORD_BIN'] = initialdir_stride
            elif sys.platform.startswith('linux') and platform.machine() == 'x86_64':
                initialdir_stride = os.path.join(os.environ['PYMOL_GIT_MOD'],"Larmord","ia64Linux2","larmord")
                os.environ['LARMORD_BIN'] = initialdir_stride
            else:
                pass
        if 'LARMORD_BIN' in os.environ:
            if VERBOSE: print 'Found LARMORD_BIN in environmental variables', os.environ['LARMORD_BIN']
            self.larmord_bin.set(os.environ['LARMORD_BIN'])
            self.larmord_para.set(os.environ['LARMORD_BIN']+"/../data/larmorD_alphas_betas_rna.dat")
            self.larmord_ref.set(os.environ['LARMORD_BIN']+"/../data/larmorD_reference_shifts_rna.dat")
        else:
            if VERBOSE: print 'LARMORD_BIN not found in environmental variables.'
            self.larmord_bin.set('')
        
        
        w = Tkinter.Label(self.dialog.interior(),text = 'LARMORD Plugin for PyMOL\nby Aaron T. Frank, 2016\n',background = 'black', foreground = 'yellow')
        w.pack(expand = 1, fill = 'both', padx = 10, pady = 5)

        # make a few tabs within the dialog
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill = 'both', expand=1, padx=10, pady=10)


        ######################
        # Tab : Structure Tab
        ######################
        page = self.notebook.add('Structure')       
        self.notebook.tab('Structure').focus_set()
        group_struc = Tkinter.LabelFrame(page, text = 'Structure')
        group_struc.pack(fill='both', expand=True, padx=10, pady=5)

        pymol_sel_ent = Pmw.EntryField(group_struc,
                                       label_text='PyMOL selection/object:',
                                       labelpos='wn',
                                       entry_textvariable=self.pymol_sel
                                       )

        larmord_bin_ent = Pmw.EntryField(group_struc,
                                      label_text='LARMORD binary:', labelpos='wn',
                                      entry_textvariable=self.larmord_bin,
                                      entry_width=20)
        larmord_bin_but = Tkinter.Button(group_struc, text = 'Browse...', command = self.getLarmordBin)


        larmord_para_ent = Pmw.EntryField(group_struc,
                                      label_text='LARMORD Parameter:', labelpos='wn',
                                      entry_textvariable=self.larmord_para,
                                      entry_width=20)
        larmord_para_but = Tkinter.Button(group_struc, text = 'Browse...', command = self.getLarmordPara)

        larmord_ref_ent = Pmw.EntryField(group_struc,
                                      label_text='LARMORD Reference Shifts:', labelpos='wn',
                                      entry_textvariable=self.larmord_ref,
                                      entry_width=20)
        larmord_ref_but = Tkinter.Button(group_struc, text = 'Browse...', command = self.getLarmordRef)

        larmord_cs_ent = Pmw.EntryField(group_struc,
                                      label_text='Chemical Shift File:', labelpos='wn',
                                      entry_textvariable=self.larmord_cs,
                                      entry_width=20)
        larmord_cs_but = Tkinter.Button(group_struc, text = 'Browse...', command = self.getLarmordCS)

        # arrange widgets using grid
        pymol_sel_ent.grid(sticky='we', row=0, column=0, columnspan=2, padx=5, pady=5)
        larmord_bin_ent.grid(sticky='we', row=3, column=0, padx=5, pady=1)
        larmord_bin_but.grid(sticky='we', row=3, column=1, padx=5, pady=1)

        larmord_para_ent.grid(sticky='we', row=4, column=0, padx=5, pady=1)
        larmord_para_but.grid(sticky='we', row=4, column=1, padx=5, pady=1)

        larmord_ref_ent.grid(sticky='we', row=5, column=0, padx=5, pady=1)
        larmord_ref_but.grid(sticky='we', row=5, column=1, padx=5, pady=1)

        larmord_cs_ent.grid(sticky='we', row=6, column=0, padx=5, pady=1)
        larmord_cs_but.grid(sticky='we', row=6, column=1, padx=5, pady=1)
        
        group_struc.columnconfigure(0, weight=9)
        group_struc.columnconfigure(1, weight=1)
      
        ######################
        # Tab : Error Tabs
        ######################
        page = self.notebook.add('Render Errors')        
        group_error = Tkinter.LabelFrame(page, text = 'Render Errors')
        group_error.grid(sticky='eswn', row=0, column=0, padx=10, pady=5)

        error_color_ent = Pmw.EntryField(group_error,
                                       label_text='Color:',
                                       labelpos='wn', value = "red",
                                       entry_textvariable=self.larmord_error_color
                                       )
        error_sel_ent = Pmw.EntryField(group_error,
                                       label_text='Type (proton or carbon?):',
                                       labelpos='wn', value = "proton",
                                       entry_textvariable=self.larmord_error_sel
                                       )
        error_color_ent.grid(sticky='we', row=0, column=0, columnspan=2, padx=5, pady=5)
        error_sel_ent.grid(sticky='we', row=1, column=0, columnspan=2, padx=5, pady=5)

        group_error.columnconfigure(0, weight=9)
        group_error.columnconfigure(1, weight=1)

        ######################
        # Tab : About Tab
        ######################
        page = self.notebook.add('About')
        group_about = Tkinter.LabelFrame(page, text = 'LARMORD Plugin for PyMOL')
        group_about.grid(sticky='we', row=0,column=0,padx=5,pady=3)
        about_plugin = """ Compute and render chemical shifts errors using LARMORD.
                   by Aaron T. Frank  <afrankz .at. umich.edu>
               Please cite this plugin if you use it in a publication.
        """
        label_about = Tkinter.Label(group_about,text=about_plugin)
        label_about.grid(sticky='we', row=0, column=0, padx=5, pady=10)        
        self.notebook.setnaturalsize()
        return

    def getLarmordBin(self):
        larmord_bin_fname = tkFileDialog.askopenfilename(title='Larmord Binary', initialdir='',filetypes=[('all','*')], parent=self.parent)
        if larmord_bin_fname: # if nonempty
            self.larmord_bin.set(larmord_bin_fname)
        return

    def getLarmordPara(self):
        larmord_para_fname = tkFileDialog.askopenfilename(title='Larmord Parameter', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_para_fname: # if nonempty
            self.larmord_para.set(larmord_para_fname)
        return

    def getLarmordRef(self):
        larmord_ref_fname = tkFileDialog.askopenfilename(title='Larmord Reference', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_ref_fname: # if nonempty
            self.larmord_ref.set(larmord_ref_fname)
        return

    def getLarmordCS(self):
        larmord_cs_fname = tkFileDialog.askopenfilename(title='Chemical Shift File', initialdir='', filetypes=[('all','*')], parent=self.parent)
        if larmord_cs_fname: # if nonempty
            self.larmord_cs.set(larmord_cs_fname)
        return

    def runLarmordOneObj(self, sel_name, objname):
        """ Run Larmord on only one object.
        
            @param one_obj_sel: the selection/object involving only one object
            @param type: string
        """
        
        one_obj_sel = '%s and %s' % (sel_name, objname)
        pdb_fn = None
        pdb_os_fh, pdb_fn = tempfile.mkstemp(suffix='.pdb') # file os handle, file name
        os.close(pdb_os_fh)
        cmd.save(filename=pdb_fn, selection=one_obj_sel)
        if VERBOSE:
            print 'Selection %s saved to %s.' % (one_obj_sel, pdb_fn)

        if pdb_fn is None:
            print 'WARNING: Larmord has no pdb file to work on!'
            return None
        
        print 'Started Running Larmord for %s ...' % (one_obj_sel,)        
        larmord_sse_dict = {}
        larmord_tmpout_os_fh, larmord_tmpout_fn = tempfile.mkstemp(suffix='.larmord')
        os.close(larmord_tmpout_os_fh)
        larmord_cmd = '%s/larmord -csfile %s -parmfile %s -reffile %s %s > %s' % (self.larmord_bin.get(), self.larmord_cs.get(), self.larmord_para.get(), self.larmord_ref.get(), pdb_fn, larmord_tmpout_fn)
        os.system(larmord_cmd)
        
        self.larmord_errors = self.parse_larmord_output(larmord_tmpout_fn)
        print 'Finished Running Larmord for %s ...' % (one_obj_sel,)
        
        return True
    
    def parse_larmord_output(self, larmord_tmpout_fn):
        """ Parse Larmord output
        
            @param larmord_tmpout_fn: larmord output file
            @param type: string
        """   
        output = {}
        parse_tmpout_os_fh, parse_tmpout_fn = tempfile.mkstemp(suffix='.larmord')        
        larmord_cmd = "awk '{print $3}' %s | tr '\n' ' ' > %s" % (larmord_tmpout_fn, parse_tmpout_fn)
        os.system(larmord_cmd)            
        fh = open(parse_tmpout_fn)
        larmord_resid = ''.join(fh.readlines()).split()

        larmord_cmd = "awk '{print $4}' %s | tr '\n' ' ' > %s" % (larmord_tmpout_fn, parse_tmpout_fn)
        os.system(larmord_cmd)            
        fh = open(parse_tmpout_fn)
        larmord_resname = ''.join(fh.readlines()).split()

        larmord_cmd = "awk '{print $5}' %s | tr '\n' ' ' > %s" % (larmord_tmpout_fn, parse_tmpout_fn)
        os.system(larmord_cmd)            
        fh = open(parse_tmpout_fn)
        larmord_nucleus = ''.join(fh.readlines()).split()

        larmord_cmd = "awk '{print $6}' %s | tr '\n' ' ' > %s" % (larmord_tmpout_fn, parse_tmpout_fn)
        os.system(larmord_cmd)            
        fh = open(parse_tmpout_fn)
        larmord_predCS = ''.join(fh.readlines()).split()

        larmord_cmd = "awk '{print $7}' %s | tr '\n' ' ' > %s" % (larmord_tmpout_fn, parse_tmpout_fn)
        os.system(larmord_cmd)            
        fh = open(parse_tmpout_fn)
        larmord_expCS = ''.join(fh.readlines()).split()
        
        cmd.alter("n. *", "b=0.0")
        
        for res in range(len(larmord_resid)):
            ch, resid, resname, nucleus, predCS, expCS = " ", larmord_resid[res], larmord_resname[res], larmord_nucleus[res], float(larmord_predCS[res]), float(larmord_expCS[res])
            k = (ch, resid, resname, nucleus)
            error = abs(predCS - expCS)
            output[k] = error
            cmd.alter("resi %s and n. %s"%(resid, nucleus), "b=%s"%error)
                 
        return output

    def render_larmord_errors(self, one_obj_sel, color="red", type = "proton" ):
        """
        """   
        cmd.bg_color( "white" )
        AllObj=cmd.get_names("all")
        cmd.set( "cartoon_ring_mode", 3 )
        cmd.set( "ray_trace_mode", 1)
        cmd.set( "ray_trace_gain", 0.1)
        cmd.set( "ray_trace_color", "black")
        cmd.set( "specular", "off")
        cmd.set( "opaque_background", "off")
        cmd.set( "antialias", 2)
        cmd.set( "all_states", 0)

        cmd.show( "cartoon")
        cmd.hide( "line")
 
        cmd.set("transparency", 0.4)
        # set VDW for all RNA atoms to 0
        cmd.alter("n. *", "vdw = 0")
        cmd.rebuild()
        # reinitialize b-factors
        stored.b = []
        # iterated over selection and store b-factors
        if (type == "proton"):
          sele = "( n. H1'+H2'+H3'+H4'+H5'+H5''+H2+H5+H6+H8)"
          stored.min = 0.0
          stored.range = 0.20   
        if (type == "carbon"):
          sele = "( n. C1'+C2'+C3'+C4'+C5'+C2+C5+C6+C8)"
          stored.min = 0.0
          stored.range = 0.90   

     
        cmd.iterate (sele, "stored.b.append(b)")
        cmd.alter(sele, "vdw = (0.4 * b) / stored.range ")
        cmd.rebuild()
        cmd.show("spheres", sele)
        cmd.set("sphere_color",color)
        cmd.set("cartoon_color","white")
 

    def runLarmord(self):
        """
        """
        # delete old results
        self.sel_obj_list = [] 
        pdb_fn = None
        sel_name= None
        sel = self.pymol_sel.get()
        
        if len(sel) > 0: # if any pymol selection/object is specified
            all_sel_names = cmd.get_names('all') # get names of all selections
            if sel in all_sel_names:
                if cmd.count_atoms(sel) == 0:
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print 'ERROR: %s' % (err_msg,)
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = sel
            # no selection/object with the input name is found
            # we assume either a single-word selector or
            # some other selection-expression is uesd
            else:
                print 'The selection/object you specified is not found.'
                print 'Your input will be interpreted as a selection-expression.'
                tmpsel = self.randomSeleName(prefix='your_sele_')
                cmd.select(tmpsel,sel)
                if cmd.count_atoms(tmpsel) == 0:
                    cmd.delete(tmpsel)
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print 'ERROR: %s' % (err_msg,)
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = tmpsel
                
        else:   # what structure do you want Larmord to work on?
            err_msg = 'No PyMOL selection/object specified!'
            print 'ERROR: %s' % (err_msg,)
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
            return False
        
        # each object in the selection is treated as an independent struc
        objlist = cmd.get_object_list(sel_name)
        self.cs_prog = 'Larmord'

        for objname in objlist:
            self.sel_obj_list.append('%s and %s' % (sel_name, objname))
            self.runLarmordOneObj(sel_name, objname)

        return True

    def runRender(self):
        """
        """
        # delete old results
        self.sel_obj_list = [] 
        pdb_fn = None
        sel_name= None
        sel = self.pymol_sel.get()
        error_color = "red"
        error_sel = "proton"
                
        if len(sel) > 0: # if any pymol selection/object is specified
            all_sel_names = cmd.get_names('all') # get names of all selections
            if sel in all_sel_names:
                if cmd.count_atoms(sel) == 0:
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print 'ERROR: %s' % (err_msg,)
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = sel
            else:
                print 'The selection/object you specified is not found.'
                print 'Your input will be interpreted as a selection-expression.'
                tmpsel = self.randomSeleName(prefix='your_sele_')
                cmd.select(tmpsel,sel)
                if cmd.count_atoms(tmpsel) == 0:
                    cmd.delete(tmpsel)
                    err_msg = 'ERROR: The selection %s is empty.' % (sel,)
                    print 'ERROR: %s' % (err_msg,)
                    tkMessageBox.showinfo(title='ERROR', message=err_msg)
                    return False
                else:
                    sel_name = tmpsel                
        else:   # what structure do you want Larmord to work on?
            err_msg = 'No PyMOL selection/object specified!'
            print 'ERROR: %s' % (err_msg,)
            tkMessageBox.showinfo(title='ERROR', message=err_msg)
            return False

        if(len(self.larmord_error_color.get()) > 0):
        		error_color = self.larmord_error_color.get()

        if(len(self.larmord_error_sel.get()) > 0):
        		error_sel = self.larmord_error_sel.get()
        
        # each object in the selection is treated as an independent struc
        objlist = cmd.get_object_list(sel_name)

        for objname in objlist:
            self.render_larmord_errors(objname, error_color, error_sel)        
        return True


    def randomSeleName(self, prefix='sele',suffix=''):
        """ generate a random selection name.
        """
        sel_list = cmd.get_names('all')
        sel_dict = dict(zip(sel_list, range(len(sel_list))))
        sel_name = '%s%d%s' % (prefix, random.randint(1000,9999), suffix)
        while(sel_name in sel_dict):
            sel_name = '%s%d%s' % (prefix, random.randint(1000,9999), suffix)
            
        return sel_name
        
        
    def execute(self, butcmd):
        """ Run the cmd represented by the botton clicked by user.
        """        
        if butcmd == 'OK':
            print 'is everything OK?'
            
        elif butcmd == 'Run LARMORD':
            rtn = self.runLarmord()
            if rtn and VERBOSE:     print 'Done with Larmord!'

        elif butcmd == 'Render Errors':
            rtn = self.runRender()
            if rtn and VERBOSE:     print 'Done rendering chemical shift errors!'

        elif butcmd == 'Exit':
            print 'Exiting LARMORD Plugin ...'
            if __name__ == '__main__':
                self.parent.destroy()
            else:
                self.dialog.withdraw()
            print 'Done.'
        else:
            print 'Exiting LARMORD Plugin because of unknown button click ...'
            self.dialog.withdraw()
            print 'Done.'
            
    
    def quit(self):
        self.dialog.destroy() 


#############################################
#
#
# Create demo in root window for testing.
#
#
##############################################
if __name__ == '__main__':
    
    class App:
        def my_show(self,*args,**kwargs):
            pass

    app = App()
    app.root = Tkinter.Tk()
    Pmw.initialise(app.root)   
    app.root.title('It seems to work!')

    widget = LARMORDPlugin(app)
    app.root.mainloop()
