package fi.tkk.ics.jbliss;

import java.util.Map;

/**
 * An interface for reporting the found generator automorphisms.
 */
public interface Reporter
{
    /**
     * The hook method that is called when a new generator automorphism
     * is found.
     *
     * @param aut         An automorphism
     * @param user_param  A parameter provided by the user
     */
    public void report(Map aut, Object user_param);
}
