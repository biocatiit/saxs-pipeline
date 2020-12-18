import platform
from ctypes import c_uint32, c_void_p, POINTER, byref
import ctypes.util
import atexit
import sys
import threading

import wx

try:
    import dbus
except Exception:
    pass

class SleepInhibit(object):
    def __init__(self):
        self.platform = platform.system()

        if self.platform == 'Darwin':
            self.sleep_inhibit = MacOSSleepInhibit()

        elif self.platform == 'Windows':
            self.sleep_inhibit = WindowsSleepInhibit()

        elif self.platform == 'Linux':
            self.sleep_inhibit = LinuxSleepInhibit()

        else:
            self.sleep_inhibit = None

        self.sleep_count = 0

    def on(self):
        if self.sleep_inhibit is not None:
            try:
                self.sleep_inhibit.on()
                self.sleep_count = self.sleep_count + 1
            except Exception:
                pass

    def off(self):
        if self.sleep_inhibit is not None:
            self.sleep_count = self.sleep_count - 1

            if self.sleep_count <= 0:
                try:
                    self.sleep_inhibit.off()
                except Exception:
                    pass

    def force_off(self):
        if self.sleep_inhibit is not None:
            try:
                self.sleep_inhibit.off()
            except Exception:
                pass

class MacOSSleepInhibit(object):
    """
    Code adapted from the python caffeine module here:
    https://github.com/jpn--/caffeine

    Used with permission under MIT license
    """

    def __init__(self):
        self.libIOKit = ctypes.cdll.LoadLibrary(ctypes.util.find_library('IOKit'))
        self.cf = ctypes.cdll.LoadLibrary(ctypes.util.find_library('CoreFoundation'))
        self.libIOKit.IOPMAssertionCreateWithName.argtypes = [ c_void_p, c_uint32, c_void_p, POINTER(c_uint32) ]
        self.libIOKit.IOPMAssertionRelease.argtypes = [ c_uint32 ]
        self.cf.CFStringCreateWithCString.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int32]
        self.cf.CFStringCreateWithCString.restype = ctypes.c_void_p

        self.kCFStringEncodingUTF8 = 0x08000100
        self._kIOPMAssertionLevelOn = 255
        self._IOPMAssertionRelease = self.libIOKit.IOPMAssertionRelease
        self._assertion = None
        self.reason = "SAXS pipeline running long process"

        self._assertID = c_uint32(0)
        self._errcode = None

        atexit.register(self.off)

    def _CFSTR(self, py_string):
        return self.cf.CFStringCreateWithCString(None, py_string.encode('utf-8'), self.kCFStringEncodingUTF8)

    def _IOPMAssertionCreateWithName(self, assert_name, assert_level, assert_msg):
        assertID = c_uint32(0)
        p_assert_name = self._CFSTR(assert_name)
        p_assert_msg = self._CFSTR(assert_msg)
        errcode = self.libIOKit.IOPMAssertionCreateWithName(p_assert_name,
            assert_level, p_assert_msg, byref(assertID))
        return (errcode, assertID)

    def _assertion_type(self, display):
        if display:
            return 'NoDisplaySleepAssertion'
        else:
            return "NoIdleSleepAssertion"

    def on(self, display=False):
        # Stop idle sleep
        a = self._assertion_type(display)
        # if a != self._assertion:
        #     self.off()
        if self._assertID.value ==0:
            self._errcode, self._assertID = self._IOPMAssertionCreateWithName(a,
        self._kIOPMAssertionLevelOn, self.reason)

    def off(self):
        self._errcode = self._IOPMAssertionRelease(self._assertID)
        self._assertID.value = 0


class WindowsSleepInhibit(object):
    """
    Prevent OS sleep/hibernate in windows; code from:
    https://github.com/h3llrais3r/Deluge-PreventSuspendPlus/blob/master/preventsuspendplus/core.py
    and
    https://trialstravails.blogspot.com/2017/03/preventing-windows-os-from-sleeping.html
    API documentation:
    https://msdn.microsoft.com/en-us/library/windows/desktop/aa373208(v=vs.85).aspx
    """

    def __init__(self):
        self.ES_CONTINUOUS = 0x80000000
        self.ES_SYSTEM_REQUIRED = 0x00000001

    def on(self):
        ctypes.windll.kernel32.SetThreadExecutionState(
            self.ES_CONTINUOUS | \
            self.ES_SYSTEM_REQUIRED)

    def off(self):
        ctypes.windll.kernel32.SetThreadExecutionState(
            self.ES_CONTINUOUS)

# For linux
class LinuxSleepInhibit(object):
    """
    Based on code from:
    https://github.com/h3llrais3r/Deluge-PreventSuspendPlus
    """
    def __init__(self):
        self.sleep_inhibitor = None
        self.get_inhibitor()

    def get_inhibitor(self):
        try:
            #Gnome session inhibitor
            self.sleep_inhibitor = GnomeSessionInhibitor()
            return
        except Exception:
            pass

        try:
            #Free desktop inhibitor
            self.sleep_inhibitor = DBusInhibitor('org.freedesktop.PowerManagement',
                '/org/freedesktop/PowerManagement/Inhibit',
                'org.freedesktop.PowerManagement.Inhibit')
            return
        except Exception:
            pass

        try:
            #Gnome inhibitor
            self.sleep_inhibitor = DBusInhibitor('org.gnome.PowerManager',
                '/org/gnome/PowerManager',
                'org.gnome.PowerManager')
            return
        except Exception:
            pass

    def on(self):
        if self.sleep_inhibitor is not None:
            self.sleep_inhibitor.inhibit()

    def off(self):
        if self.sleep_inhibitor is not None:
            self.sleep_inhibitor.uninhibit()

class DBusInhibitor:
    def __init__(self, name, path, interface, method=['Inhibit', 'UnInhibit']):
        self.name = name
        self.path = path
        self.interface_name = interface

        bus = dbus.SessionBus()
        devobj = bus.get_object(self.name, self.path)
        self.iface = dbus.Interface(devobj, self.interface_name)
        # Check we have the right attributes
        self._inhibit = getattr(self.iface, method[0])
        self._uninhibit = getattr(self.iface, method[1])

    def inhibit(self):
        self.cookie = self._inhibit('SAXS Pipeline', 'long_process')

    def uninhibit(self):
        self._uninhibit(self.cookie)


class GnomeSessionInhibitor(DBusInhibitor):
    TOPLEVEL_XID = 0
    INHIBIT_SUSPEND = 4

    def __init__(self):
        DBusInhibitor.__init__(self, 'org.gnome.SessionManager',
                               '/org/gnome/SessionManager',
                               'org.gnome.SessionManager',
                               ['Inhibit', 'Uninhibit'])

    def inhibit(self):
        self.cookie = self._inhibit('SAXS Pipeline',
                                    GnomeSessionInhibitor.TOPLEVEL_XID,
                                    'long_process',
                                    GnomeSessionInhibitor.INHIBIT_SUSPEND)


def signal_handler(sig, frame):
    main_frame = wx.Window.FindWindowByName('PipelineFrame')
    main_frame.cleanup_and_quit_forced()


class MyModalDialogHook(wx.ModalDialogHook):

    def __init__(self):
        wx.ModalDialogHook.__init__(self)

        self.all_modal_dialogs = []

    def Enter(self, dialog):
        self.all_modal_dialogs.append(dialog)

        return wx.ID_NONE

    def Exit(self, dialog):
        if dialog in self.all_modal_dialogs:
            self.all_modal_dialogs.remove(dialog)


def setup_thread_excepthook():
    """
    Workaround for `sys.excepthook` thread bug from:
    http://bugs.python.org/issue1230540

    Call once from the main thread before creating any threads.
    """

    init_original = threading.Thread.__init__

    def init(self, *args, **kwargs):

        init_original(self, *args, **kwargs)
        run_original = self.run

        def run_with_except_hook(*args2, **kwargs2):
            try:
                run_original(*args2, **kwargs2)
            except Exception:
                sys.excepthook(*sys.exc_info())

        self.run = run_with_except_hook

    threading.Thread.__init__ = init
