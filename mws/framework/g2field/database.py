# Code for interactive with google sheets that contain shim data.
import httplib2
import os

from apiclient import discovery
import oauth2client
from oauth2client import client
from oauth2client import tools

import gspread
import time

flags = None

# If modifying these scopes, delete your previously saved credentials
# at ~/.credentials/drive-python-quickstart.json
SCOPES = 'https://spreadsheets.google.com/feeds'
CLIENT_SECRET_FILE = 'client_secret.json'
APPLICATION_NAME = 'Tilt Data Cruncher'
SEC_PER_DAY = 3600 * 24
DATA_DIR = 'data/'
START_DATE = "01-02-2016"
TOP_HAT_SETTINGS_BOOK = 'Settings - Top Hats'
TOP_HAT_SETTINGS_SHEET = 'Current'
WEDGE_SETTINGS_BOOK = 'Settings - Wedges'
WEDGE_SETTINGS_SHEET = 'Current'

wedge_turns_per_mm = 1.6
N_TOP_HATS = 48
N_WEDGES = 864


def get_gspread_session():
    """Get an authorized handle to grab data from the google sheets."""
    credentials = get_credentials()
    return gspread.authorize(credentials)


def get_credentials():
    """Gets valid user credentials from storage.

    If nothing has been stored, or if the stored credentials are invalid,
    the OAuth2 flow is completed to obtain the new credentials.

    Returns:
        Credentials, the obtained credential.
    """
    home_dir = os.path.expanduser('~')
    credential_dir = os.path.join(home_dir, '.credentials')
    if not os.path.exists(credential_dir):
        os.makedirs(credential_dir)
    credential_path = os.path.join(credential_dir,
                                   'drive-python-tilt-data.json')

    store = oauth2client.file.Storage(credential_path)
    credentials = store.get()
    if not credentials or credentials.invalid:
        flow = client.flow_from_clientsecrets(CLIENT_SECRET_FILE, SCOPES)
        flow.user_agent = APPLICATION_NAME
        if flags:
            credentials = tools.run_flow(flow, store, flags)
        else: # Needed only for compatibility with Python 2.6
            credentials = tools.run(flow, store)
        print('Storing credentials to ' + credential_path)
    return credentials


def upload_top_hat_changes(deltas, datafile, units='mm'):
    """Upload the current adjustment plan to the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    book = gc.open(TOP_HAT_SETTINGS_BOOK)
    sheet = book.worksheet(TOP_HAT_SETTINGS_SHEET)

    cells = sheet.range('A1:G49')

    if len(deltas) == N_TOP_HATS / 2:

        deltas = np.hstack((deltas, deltas))

    if units == 'mm':

        deltas /= 0.0254

    elif units == 'mils':

        pass

    else:
        raise ValueError("Not a recognized unit.")

    cells[0].value = 'Yoke ID'
    cells[1].value = 'Shim ID'
    cells[2].value = 'Layer'
    cells[3].value = 'Expected'
    cells[4].value = 'Change'
    cells[5].value = 'Input'
    cells[6].value = 'Modified'

    for i in xrange(1, N_TOP_HATS + 1):

        cells[7 * i + 0].value = chr(ord('A') + ((i - 1) % 24) / 2)
        cells[7 * i + 1].value = 1 + (i - 1) % 2

        if i > N_TOP_HATS / 2:
            cells[7 * i + 2].value = 'BOT'
        else:
            cells[7 * i + 2].value = 'TOP'

        cells[7 * i + 4].value = '{:.2f}'.format(deltas[i - 1])
        cells[7 * i + 5].value = os.path.basename(datafile)
        cells[7 * i + 6].value = time.strftime('%m/%d/%Y')

    sheet.update_cells(cells)


def upload_wedge_changes(deltas, datafile, units='mm'):
    """Upload the current adjustment plan to the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    book = gc.open(WEDGE_SETTINGS_BOOK)
    sheet = book.worksheet(WEDGE_SETTINGS_SHEET)

    cells = sheet.range('A1:G865')

    if len(deltas) == N_WEDGES / 2:

        d = np.empty(N_WEDGES)
        for i in xrange(N_WEDGES):
            d[i] = deltas[i / 2]

        deltas = np.array(d)

    if units == 'mm':

        deltas /= 1.6

    elif units == 'turns':

        pass

    else:
        raise ValueError("Not a recognized unit.")

    cells[0].value = 'Pole ID'
    cells[1].value = 'Shim ID'
    cells[2].value = 'Layer'
    cells[3].value = 'Expected'
    cells[4].value = 'Change'
    cells[5].value = 'Input'
    cells[6].value = 'Modified' 

    for i in xrange(1, N_WEDGES + 1):

        cells[7 * i + 0].value = (i - 1) / 24 + 1
        cells[7 * i + 1].value = 1 + ((i - 1) / 2) % 12

        if i % 2 == 1:
            cells[7 * i + 2].value = 'TOP'
        else:
            cells[7 * i + 2].value = 'BOT'

        cells[7 * i + 4].value = '{:.2f}'.format(deltas[i - 1])
        cells[7 * i + 5].value = os.path.basename(datafile)
        cells[7 * i + 6].value = time.strftime('%m/%d/%Y')

    sheet.update_cells(cells)


def get_current_top_hat_cells(sheetname):
    """Loads the cells the current top hat google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    try:
        top_hat_sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(sheetname)

    except(gspread.WorksheetNotFound):
        print "Could not find 'Current Top Hats' worksheet.  Creating blank one."
        sh = gc.open(SHIM_SETTINGS_SHEET)
        top_hat_sheet = sh.add_worksheet(sheetname, rows="13", cols="6")
    
    return top_hat_sheet.range('A1:F13')


def get_current_top_hat_settings(sheetname=TOP_HAT_SETTINGS_SHEET):
    """Loads and returns the current top hat settings as an array."""
    top_hat_cells = get_current_top_hat_cells(sheetname)

    top_hat_pos = np.zeros([12, 2])
    
    for i in xrange(12):
        for j in xrange(2):
            top_hat_pos[i, j] += float(top_hat_cells[6 * i + 2 * j + 7].value)
            top_hat_pos[i, j] += float(top_hat_cells[6 * i + 2 * j + 8].value)
            top_hat_pos[i, j] *= 0.5
    
    return top_hat_pos.reshape(24) * 0.0254


def update_top_hat_settings(top_hat_deltas, sheetname=TOP_HAT_SETTINGS_SHEET):
    """Update the current positions of the top hats in the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)

    # Grab the sheet for updating top hats, create if need be.
    try:
        top_hat_sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(NEW_TOP_HAT_SETTINGS_SHEET)
        
    except(gspread.WorksheetNotFound):
        sh = gc.open(SHIM_SETTINGS_SHEET)
        top_hat_sheet = sh.add_worksheet(title=NEW_TOP_HAT_SETTINGS_SHEET, rows="13", cols="6")

    top_hat_cells = get_current_top_hat_cells(sheetname)
    new_cells = top_hat_sheet.range('A1:F13')

    delta = top_hat_deltas.reshape(12, 2)

    # Set the column headers.
    top_hat_cells[0].value = "Yoke"
    top_hat_cells[1].value = "BOT-1"
    top_hat_cells[2].value = "TOP-1"
    top_hat_cells[3].value = "BOT-2"
    top_hat_cells[4].value = "TOP-2"
    top_hat_cells[5].value = "Last Modified"
    
    # Fill the values.
    for i in xrange(12):
        for j in xrange(2):
            val = float(top_hat_cells[6 * i + 2 * j + 7].value)
            top_hat_cells[6 * i + 2 * j + 7].value = str(val + delta[i, j] / 0.0254)

            val = float(top_hat_cells[6 * i + 2 * j + 8].value)
            top_hat_cells[6 * i + 2 * j + 8].value = str(val + delta[i, j] / 0.0254)
        
        ws_date = time.strftime('%m/%d/%Y', time.localtime(int(time.time())))
        top_hat_cells[6 * i + 11].value = ws_date
    
    for i in xrange(len(top_hat_cells)):
        new_cells[i].value = top_hat_cells[i].value
        
    top_hat_sheet.update_cells(new_cells)
    
    
def get_current_wedge_cells(sheetname):
    """Loads the cells the current wegdge shim google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)
    
    try:
        sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(sheetname)

    except(gspread.WorksheetNotFound):
        print "Could not find 'Current Wedges' worksheet.  Creating blank one."
        sh = gc.open(SHIM_SETTINGS_SHEET)
        sheet = sh.add_worksheet(sheetname, rows="73", cols="14")
    
    return sheet.range('A1:N73')


def get_current_wedge_settings(sheetname=WEDGE_SETTINGS_SHEET):
    """Loads and returns the current top hat settings as an array."""
    wedge_cells = get_current_wedge_cells(sheetname)

    wedge_pos = np.zeros([72, 12])
    
    for i in xrange(72):
        for j in xrange(12):
            val = wedge_cells[14 * i + j + 15].value
            
            if val == '':
                wedge_pos[i, j] = 0.0
                
            else:
                wedge_pos[i, j] = float(wedge_cells[14 * i + j + 15].value)

    return wedge_pos.reshape(72 * 12) / wedge_turns_per_mm

    
def update_wedge_settings(wedge_deltas, sheetname=WEDGE_SETTINGS_SHEET):
    """Update the current positions of the top hats in the google sheet."""
    credentials = get_credentials()
    gc = gspread.authorize(credentials)

    # Grab the sheet for updating top hats, create if need be.
    try:
        wedge_sheet = gc.open(SHIM_SETTINGS_SHEET).worksheet(NEW_WEDGE_SETTINGS_SHEET)
        
    except(gspread.WorksheetNotFound):
        sh = gc.open(SHIM_SETTINGS_SHEET)
        wedge_sheet = sh.add_worksheet(title=NEW_WEDGE_SETTINGS_SHEET, rows="73", cols="14")

    wedge_cells = get_current_wedge_cells(sheetname)
    new_cells = wedge_sheet.range('A1:N73')

    delta = wedge_deltas.reshape([72, 12])
    
    # Set the column headers.
    wedge_cells[0].value = "Pole ID"
    
    for i in xrange(1, 13):
        wedge_cells[i].value = "W%02i" % i

    wedge_cells[13].value = "Last Modified"
    
    # Set the values.
    for i in xrange(72):
        for j in xrange(12):
            idx = 14 * i + j + 15

            val = wedge_cells[idx].value
            if val == '': 
                val = 0.0
            else:
                val = float(val)
            
            wedge_cells[idx].value = str(val + delta[i, j] / wedge_turns_per_mm)

        # Set the Pole ID.
        if i / 36 == 0:
            wedge_cells[14 * i + 14].value = 'BOT-%03i' % (i % 36 + 1)

        else:
            wedge_cells[14 * i + 14].value = 'TOP-%03i' % (i % 36 + 1)

        # Set the modification date.
        ws_date = time.strftime('%m/%d/%Y', time.localtime(int(time.time())))
        wedge_cells[14 * i + 27].value = ws_date
    
    for i in xrange(len(wedge_cells)):
        new_cells[i].value = wedge_cells[i].value
        
    wedge_sheet.update_cells(new_cells)

